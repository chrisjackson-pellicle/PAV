#!/usr/bin/env python3
"""
Build per-order plastid genome collections for angiosperms.
- Reads angiosperm orders from plastid_annotation_validator/data/all_orders.txt
- Selects up to 15 unique accessions per order (favoring family diversity)
- Writes genomes and metadata to plastid_annotation_validator/data/order_genomes/<Order>/
- Requires Entrez Direct tools: esearch, efetch, esummary, xtract

Notes:
- Uses a single nucleotide DB query to find plastid/chloroplast complete genomes; no RefSeq-only restriction
- Attempts to maximize unique families, then fills remaining slots
- If fewer than 15 are available for an order, uses as many as exist
"""

import os
import sys
import shlex
import time
import json
import subprocess
from collections import defaultdict, OrderedDict
from pathlib import Path
from typing import Dict, List, Tuple, Set

REQUIRED_CMDS = ["esearch", "efetch", "esummary", "xtract"]

# Global politeness delay between downloads (seconds)
DOWNLOAD_DELAY_SECONDS = 0.6
# Max records to pull per order to select 15 diverse entries
RETMAX_PER_ORDER = 400


def check_dependencies() -> None:
	missing = []
	for cmd in REQUIRED_CMDS:
		try:
			result = subprocess.run([cmd, "-help"], capture_output=True, text=True)
			if result.returncode != 0:
				missing.append(cmd)
		except FileNotFoundError:
			missing.append(cmd)
	if missing:
		sys.stderr.write(
			"Missing required Entrez Direct tools: {}\nInstall via conda: conda install -c bioconda entrez-direct\n".format(
				", ".join(missing)
			)
		)
		sys.exit(1)


def run_cmd(cmd: List[str]) -> str:
	proc = subprocess.run(cmd, capture_output=True, text=True)
	if proc.returncode != 0:
		raise RuntimeError(f"Command failed: {' '.join(shlex.quote(c) for c in cmd)}\nSTDERR: {proc.stderr.strip()}")
	return proc.stdout


def run_pipeline(command_str: str) -> str:
	"""Run a shell pipeline and return stdout text (raises on non-zero)."""
	proc = subprocess.run(["bash", "-lc", command_str], capture_output=True, text=True)
	if proc.returncode != 0:
		raise RuntimeError(f"Pipeline failed: {command_str}\nSTDERR: {proc.stderr.strip()}")
	return proc.stdout


def discover_angiosperm_orders() -> List[str]:
	"""Discover order-rank taxa under flowering plants (Magnoliophyta) via taxonomy DB."""
	# Magnoliophyta taxid (angiosperms) is 3398; use subtree and filter by Rank==order
	pipeline = (
		"esearch -db taxonomy -query 'txid3398[Subtree]' | "
		"efetch -db taxonomy -format xml | "
		"xtract -pattern Taxon -if Rank -equals order -element ScientificName"
	)
	out = run_pipeline(pipeline)
	orders = sorted({ln.strip() for ln in out.splitlines() if ln.strip()})
	return orders


def read_orders_from_file(orders_file: Path) -> List[str]:
	"""Read order names from a plain text file (one per line, '#' comments allowed)."""
	with open(orders_file, "r") as fh:
		orders = [ln.strip() for ln in fh if ln.strip() and not ln.strip().startswith("#")]
	return orders


def get_order_taxid(order_name: str) -> str:
	"""Get taxonomy TaxId for a given order name.

	Strategy:
	- Search taxonomy by scientific name scoped to angiosperms subtree (txid3398)
	- Filter by Rank==order during xtract
	- Fallback: retry without subtree constraint
	"""
	# Primary: scoped to angiosperms
	query1 = f'"{order_name}"[SCIN] AND txid3398[Subtree]'
	pipeline1 = (
		f"esearch -db taxonomy -query {shlex.quote(query1)} | "
		"efetch -db taxonomy -format xml | "
		"xtract -pattern Taxon -if Rank -equals order -element TaxId"
	)
	try:
		out1 = run_pipeline(pipeline1).strip()
		if out1:
			return out1.splitlines()[0].strip()
	except RuntimeError:
		pass

	# Fallback: name only, still enforce Rank==order in xtract
	query2 = f'"{order_name}"[SCIN]'
	pipeline2 = (
		f"esearch -db taxonomy -query {shlex.quote(query2)} | "
		"efetch -db taxonomy -format xml | "
		"xtract -pattern Taxon -if Rank -equals order -element TaxId"
	)
	try:
		out2 = run_pipeline(pipeline2).strip()
		if out2:
			return out2.splitlines()[0].strip()
	except RuntimeError:
		pass

	return ""


def fetch_order_candidates(order_name: str) -> List[Tuple[str, str, str]]:
	"""Return list of (accession, organism, taxid) candidates for an order.
	Use a single nucleotide DB query focusing on plastid/chloroplast complete genomes.
	"""
	def query_candidates(q: str) -> List[Tuple[str, str, str]]:
		pipe = (
			f"esearch -db nucleotide -query {shlex.quote(q)} | "
			"efetch -db nucleotide -format docsum | "
			"xtract -pattern DocumentSummary -element AccessionVersion Organism TaxId"
		)
		# print(f"  pipe: {pipe}")
		try:
			text = run_pipeline(pipe)
			# print(f"  text: {text}")
		except RuntimeError:
			return []
		out: List[Tuple[str, str, str]] = []
		for line in text.splitlines():
			parts = [p.strip() for p in line.split("\t")]
			if len(parts) == 3 and all(parts):
				out.append((parts[0], parts[1], parts[2]))
		return out

	taxid = get_order_taxid(order_name)

	print(f"Order taxid: {taxid}")
	if not taxid:
		return []
	# Single query: nucleotide DB, plastid/chloroplast filters, complete genome, organism expanded to subtree taxa
	query = (
		f"txid{taxid}[Organism:exp] AND (chloroplast[Filter] OR plastid[Filter] OR plastome[Filter]) AND \"complete genome\""
	)
	candidates = query_candidates(query)
	print(f"  nucleotide_single_query count: {len(candidates)}")
	return candidates


def fetch_tax_ranks_for_taxids(taxids: List[str]) -> Dict[str, Dict[str, str]]:
	"""Return mapping taxid -> {rank_name: scientific_name} for ranks in LineageEx (including family, order)."""
	result: Dict[str, Dict[str, str]] = {}
	unique_taxids = list(OrderedDict.fromkeys(taxids))
	print(f"  unique_taxids: {unique_taxids}")
	chunk_size = 200
	for i in range(0, len(unique_taxids), chunk_size):
		chunk = unique_taxids[i : i + chunk_size]
		if not chunk:
			continue
		ids_csv = ",".join(chunk)
		# Map TaxId -> family/order using efetch XML for reliable LineageEx parsing
		fam_txt = run_pipeline(
			f"efetch -db taxonomy -id {shlex.quote(ids_csv)} -format xml | "
			"xtract -pattern Taxon -element TaxId -block LineageEx -if Rank -equals family -element ScientificName"
		)
		ord_txt = run_pipeline(
			f"efetch -db taxonomy -id {shlex.quote(ids_csv)} -format xml | "
			"xtract -pattern Taxon -element TaxId -block LineageEx -if Rank -equals order -element ScientificName"
		)
		# print(f"  fam_txt: {fam_txt}")
		# print(f"  ord_txt: {ord_txt}")
		def extract_family_from_parts(parts: List[str]) -> str:
			for token in parts[1:]:
				if token.lower().endswith("aceae"):
					return token
			return ""

		def extract_order_from_parts(parts: List[str]) -> str:
			for token in parts[1:]:
				if token.lower().endswith("ales"):
					return token
			return ""

		for line in fam_txt.splitlines():
			parts = [p.strip() for p in line.split("\t") if p.strip()]
			if not parts:
				continue
			taxid_key = parts[0]
			if not taxid_key.isdigit():
				continue
			family_name = extract_family_from_parts(parts)
			if family_name:
				result.setdefault(taxid_key, {})["family"] = family_name

		for line in ord_txt.splitlines():
			parts = [p.strip() for p in line.split("\t") if p.strip()]
			if not parts:
				continue
			taxid_key = parts[0]
			if not taxid_key.isdigit():
				continue
			order_name = extract_order_from_parts(parts)
			if order_name:
				result.setdefault(taxid_key, {})["order"] = order_name
	
	# print(f"  result: {result}")
	return result


def select_diverse_by_family(candidates: List[Tuple[str, str, str]], tax_ranks: Dict[str, Dict[str, str]], max_per_order: int = 15) -> List[Tuple[str, str, str, str]]:
	"""
	Pick up to max_per_order entries maximizing unique families first.
	Avoid repeating species across accessions.
	Return list of (accession, organism, family, order).
	"""

	def normalize_species_for_dedup(organism_name: str) -> str:
		name = organism_name.strip().lower()
		parts = [p for p in name.split() if p]
		if len(parts) >= 2:
			return f"{parts[0]} {parts[1]}"
		return name

	selected: List[Tuple[str, str, str, str]] = []
	seen_accessions: Set[str] = set()
	used_families: Set[str] = set()
	seen_species: Set[str] = set()
	order_name_from_tax: str = ""
	# First pass: unique families
	for acc, org, taxid in candidates:
		ranks = tax_ranks.get(taxid, {})
		family = ranks.get("family", "")
		order_name = ranks.get("order", "")
		if not order_name_from_tax and order_name:
			order_name_from_tax = order_name
		species_key = normalize_species_for_dedup(org)
		if family and acc not in seen_accessions and family not in used_families and species_key not in seen_species:
			selected.append((acc, org, family, order_name))
			seen_accessions.add(acc)
			used_families.add(family)
			seen_species.add(species_key)
			if len(selected) >= max_per_order:
				return selected
	# Second pass: fill remaining even if families repeat (but still avoid species repeats)
	for acc, org, taxid in candidates:
		if len(selected) >= max_per_order:
			break
		if acc in seen_accessions:
			continue
		ranks = tax_ranks.get(taxid, {})
		family = ranks.get("family", "")
		order_name = ranks.get("order", "")
		species_key = normalize_species_for_dedup(org)
		if species_key in seen_species:
			continue
		selected.append((acc, org, family, order_name))
		seen_accessions.add(acc)
		seen_species.add(species_key)
	return selected


def write_order_metadata(order_dir: Path, order_name: str, selections: List[Tuple[str, str, str, str]]) -> None:
	# TSV
	tsv_path = order_dir / "plastid_genomes_taxonomy.tsv"
	with open(tsv_path, "w") as f:
		f.write("Accession\tSpecies\tCommon_Name\tFamily\tOrder\tFile_Path\n")
		for acc, organism, family, ord_name in selections:
			# File path is relative to data/ directory for portability
			f.write(f"{acc}\t{organism}\t\t{family}\t{order_name}\torder_genomes/{order_dir.name}/{acc}.gb\n")
	# README
	readme_path = order_dir / "README.md"
	family_to_species: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
	for acc, organism, family, _ in selections:
		family_to_species[family or "Unknown"].append((organism, acc))
	lines: List[str] = []
	lines.append(f"# {order_name} Plastid Reference Genomes\n\n")
	lines.append("This directory contains plastid (chloroplast) genomes from the order, selected to maximize family diversity.\n\n")
	lines.append("## Collection Summary\n\n")
	lines.append(f"- Total Genomes: {len(selections)}\n")
	lines.append(f"- Order: {order_name}\n")
	lines.append(f"- Families Represented: {len([k for k in family_to_species.keys() if k and k != 'Unknown'])}\n\n")
	lines.append("## Families and Species\n\n")
	for fam in sorted(family_to_species.keys()):
		lines.append(f"### {fam or 'Unknown'}\n")
		for organism, acc in sorted(family_to_species[fam]):
			lines.append(f"- {organism} — {acc}\n")
		lines.append("\n")
	with open(readme_path, "w") as f:
		f.write("".join(lines))


def download_accessions(order_dir: Path, selections: List[Tuple[str, str, str, str]]) -> Tuple[int, int]:
	success, fail = 0, 0
	for acc, _, _, _ in selections:
		out_file = order_dir / f"{acc}.gb"
		if out_file.exists():
			success += 1
			continue
		cmd = [
			"efetch", "-db", "nucleotide", "-id", acc, "-format", "gb"
		]
		try:
			with open(out_file, "w") as fh:
				proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
			if proc.returncode == 0 and out_file.exists() and out_file.stat().st_size > 0:
				success += 1
			else:
				fail += 1
		except Exception:
			fail += 1
		time.sleep(DOWNLOAD_DELAY_SECONDS)
	return success, fail


def main() -> None:
	check_dependencies()
	# Resolve data locations relative to this script
	script_dir = Path(__file__).resolve().parent
	# Determine the 'data' directory whether this script is in repo root or inside data/download_scripts
	candidate_data_dirs = [
		script_dir.parent,  # when script is inside .../data/download_scripts
		script_dir / "plastid_annotation_validator" / "data",  # when script is at repo root
	]
	data_dir: Path
	for cand in candidate_data_dirs:
		if (cand / "all_orders.txt").exists():
			data_dir = cand
			break
	else:
		# Fallback: assume parent of script dir is the data directory
		data_dir = script_dir.parent

	orders_file = data_dir / "all_orders.txt"
	output_root = data_dir / "order_genomes"
	output_root.mkdir(parents=True, exist_ok=True)

	print(f"Reading angiosperm orders from file: {orders_file}")
	orders = read_orders_from_file(orders_file)
	print(f"Found {len(orders)} orders to process")

	for idx, order_name in enumerate(orders, start=1):
		print(f"\n[{idx}/{len(orders)}] Processing order: {order_name}")
		try:
			candidates = fetch_order_candidates(order_name)
			if not candidates:
				print(f"  No plastid chloroplast complete genomes found for {order_name}")
				continue
			# Build taxonomy ranks map for all candidate taxids
			taxids = [t for _, _, t in candidates]
			print(f"  taxids: {taxids}")
			ranks = fetch_tax_ranks_for_taxids(taxids)
			selections = select_diverse_by_family(candidates, ranks, max_per_order=15)
			if not selections:
				print(f"  No selections possible for {order_name}")
				continue
			# Create directory and metadata under data/order_genomes/<Order>
			order_dir = output_root / order_name
			order_dir.mkdir(exist_ok=True)
			write_order_metadata(order_dir, order_name, selections)
			# Download
			print(f"  Downloading {len(selections)} accessions for {order_name}...")
			succ, fail = download_accessions(order_dir, selections)
			print(f"  Completed {order_name}: success={succ}, failed={fail}")
		except Exception as e:
			print(f"  Error processing {order_name}: {e}")

	print("\nAll done.")


if __name__ == "__main__":
	main()


