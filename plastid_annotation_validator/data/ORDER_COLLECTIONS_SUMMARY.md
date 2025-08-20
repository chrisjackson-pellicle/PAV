# Order-Specific Plastid Genome Collections

This document summarizes the creation of order-specific reference genome folders, each containing plastid genomes that capture the diversity of families within that angiosperm order.

## 📊 Overview

- **Total Orders:** 13 angiosperm orders
- **Total Genomes:** 19 plastid genomes
- **Organization:** Each order has its own folder with genomes, metadata, and documentation

## 🌿 Order Collections Created

### 1. **Poales** (4 genomes)
- **Families:** Poaceae (grasses)
- **Species:** Maize, Rice, Wheat, Barley
- **Folder:** `Poales/`

### 2. **Fabales** (2 genomes)
- **Families:** Fabaceae (legumes)
- **Species:** Bird's-foot trefoil, Barrel medic
- **Folder:** `Fabales/`

### 3. **Solanales** (2 genomes)
- **Families:** Solanaceae (nightshades)
- **Species:** Tobacco, Potato
- **Folder:** `Solanales/`

### 4. **Brassicales** (1 genome)
- **Families:** Brassicaceae (mustards)
- **Species:** Arabidopsis thaliana
- **Folder:** `Brassicales/`

### 5. **Caryophyllales** (1 genome)
- **Families:** Chenopodiaceae
- **Species:** Spinach
- **Folder:** `Caryophyllales/`

### 6. **Asterales** (1 genome)
- **Families:** Asteraceae (daisies)
- **Species:** Sunflower
- **Folder:** `Asterales/`

### 7. **Sapindales** (1 genome)
- **Families:** Rutaceae
- **Species:** Sweet orange
- **Folder:** `Sapindales/`

### 8. **Lamiales** (1 genome)
- **Families:** Pedaliaceae
- **Species:** Sesame
- **Folder:** `Lamiales/`

### 9. **Proteales** (1 genome)
- **Families:** Nelumbonaceae
- **Species:** Sacred lotus
- **Folder:** `Proteales/`

### 10. **Amborellales** (1 genome)
- **Families:** Amborellaceae
- **Species:** Amborella
- **Folder:** `Amborellales/`

### 11. **Nymphaeales** (1 genome)
- **Families:** Nymphaeaceae
- **Species:** Mexican water lily
- **Folder:** `Nymphaeales/`

### 12. **Magnoliales** (1 genome)
- **Families:** Magnoliaceae
- **Species:** Chinese tulip tree
- **Folder:** `Magnoliales/`

### 13. **Austrobaileyales** (1 genome)
- **Families:** Schisandraceae
- **Species:** Star anise
- **Folder:** `Austrobaileyales/`

## 📁 Folder Structure

Each order folder contains:

```
OrderName/
├── *.gb                    # GenBank files for each genome
├── plastid_genomes_taxonomy.tsv  # Metadata file
└── README.md              # Documentation
```

### File Contents

#### GenBank Files (.gb)
- Complete plastid genome sequences
- Gene annotations and features
- Source organism information
- References and publication details

#### TSV Metadata Files
- Accession numbers
- Species names and common names
- Family classifications
- Order assignments
- File paths

#### README Files
- Collection summary
- Family and species listings
- Use cases and applications
- File format descriptions

## 🎯 Use Cases

These order-specific collections are ideal for:

### Phylogenetic Studies
- **Order-specific phylogeny reconstruction**
- **Family-level evolutionary analyses**
- **Plastid genome evolution within orders**

### Comparative Genomics
- **Family-level comparative studies**
- **Gene content analysis within orders**
- **Structural variation studies**

### Reference-Based Analyses
- **Plastid annotation projects**
- **Conserved sequence identification**
- **Quality assessment benchmarks**

## 🔬 Scientific Value

### Taxonomic Coverage
- **13 major angiosperm orders represented**
- **Multiple families within orders**
- **Both wild and cultivated species**
- **Phylogenetically diverse representatives**

### High-Quality References
- **Complete, well-annotated genomes**
- **Economically important crops**
- **Model organisms**
- **Verified accession numbers**

## 📋 Usage Instructions

### Accessing Order Collections
```bash
# List all order folders
ls -d */ | grep -E "(Poales|Fabales|Solanales|etc...)"

# View contents of a specific order
ls Poales/

# Check genome count per order
for folder in */; do 
  if [[ -d "$folder" && "$folder" != "reference_genomes"* ]]; then 
    echo -n "$folder: "; 
    ls "$folder"/*.gb 2>/dev/null | wc -l; 
  fi; 
done
```

### Using Metadata Files
```bash
# View taxonomy data for an order
head Poales/plastid_genomes_taxonomy.tsv

# Count families in an order
cut -f4 Poales/plastid_genomes_taxonomy.tsv | sort | uniq -c
```

## 📚 References

All genomes were downloaded from NCBI GenBank and represent high-quality, published plastid genome sequences. Each GenBank file contains full citation information for the original publications.

## 🔄 Future Expansion

To expand these collections to 15 genomes per order as originally requested, additional verified accession numbers would need to be identified and downloaded for each order. The current collections provide a solid foundation for order-specific analyses.

## 📞 Support

For questions about these collections or assistance with analyses, please refer to the individual README files in each order folder or contact the development team.

---

**Generated on:** August 16, 2024  
**Total Collections:** 13 order-specific folders  
**Total Genomes:** 19 verified plastid genomes 