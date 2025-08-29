# A Bioinformatics Approach to Studying Convergent Evolution in Cancer

## Overview
This master's project develops a computational model based on the Moran process to simulate tumor evolution and compare two fundamental paradigms in cancer biology: **oncogene addiction** versus **convergent evolution**. The study challenges traditional "one gene, one drug" therapeutic approaches by demonstrating how multiple genetic pathways can converge on identical survival advantages.

## Research Question
Can computational modeling differentiate between oncogene addiction and convergent evolution paradigms in cancer, and what are the implications for therapeutic strategy design?

## Key Hypothesis
- **Oncogene addiction**: Genetic and phenotypic diversity should decline together (coupled evolution)
- **Convergent evolution**: Genetic diversity should be maintained while phenotypic diversity decreases (decoupled evolution)

## Biological Context
Traditional cancer treatment assumes **oncogene addiction** - that specific genes drive cancer through deterministic pathways, leading to targeted therapies. However, clinical evidence suggests cancer follows **evolutionary principles** where:
- Multiple genetic routes can achieve identical phenotypic outcomes (drug resistance, proliferation)
- Tumor populations maintain genetic heterogeneity while converging on survival advantages
- Single-target therapies often fail due to evolutionary escape through alternative pathways

## Technologies Used
- **Python** - Primary simulation language
- **NumPy** - Numerical computations and population modeling
- **Matplotlib/Seaborn** - Data visualization
- **SciPy** - Statistical analysis (Shannon entropy, Pearson correlation)
- **Moran Process** - Stochastic population genetics model
- **Shannon Entropy** - Quantifying genetic and phenotypic diversity

## Computational Model Design

### Population Structure
- **Population size**: 200 cells
- **Genome representations**: 3-bit (8 genotypes) and 12-bit (4,096 genotypes)
- **Fitness calculation**: Additive model (sum of 1s in binary genotype)
- **Selection pressure**: Deterministic replacement of least fit with most fit
- **Mutation rate**: 0.01 per bit per generation

### Model Variations
1. **3-bit genome, 50 generations** - Limited genotypic space
2. **12-bit genome, 50 generations** - Expanded evolutionary possibilities  
3. **12-bit genome, 200 generations** - Extended temporal dynamics
4. **Peak fitness model** - Convergent evolution with intermediate fitness optimum

## Key Findings

### Quantitative Results
| Model Type | Genotypic-Phenotypic Correlation | Final Genetic Diversity | Final Phenotypic Diversity |
|------------|--------------------------------|------------------------|---------------------------|
| 3-bit Oncogene | Strong coupling (~1.6 ratio) | 2.6 bits | 1.6 bits |
| 12-bit Oncogene | Moderate coupling | 6.5 bits | 2.7 bits |
| Convergent Evolution | **Decoupled** | **High maintenance** | **Reduced variation** |

### Critical Discovery
**12-bit genome simulations demonstrated clear decoupling of genetic and phenotypic diversity** - genotypic diversity declined 17% while phenotypic diversity remained stable (7% decline), providing computational evidence for convergent evolution in cancer.

## Clinical Implications

### Challenge to Current Paradigms
- **Oncogene addiction limitations**: Single-target therapies may be inherently limited by evolutionary flexibility
- **Treatment resistance mechanisms**: Multiple genetic pathways can achieve identical resistance phenotypes
- **Therapeutic strategy shift**: From "one gene, one drug" to evolution-guided combination approaches

### Evidence-Based Recommendations
1. **Combination therapies** targeting multiple evolutionary pathways simultaneously
2. **Adaptive treatment strategies** that account for tumor evolutionary potential
3. **Kinetic rather than genetic** biomarkers for treatment selection
4. **Evolutionary modeling** for predicting treatment outcomes (clinical R² = 0.81)

## Files in this Repository

### Analysis Code
- `Convergent_Evolution.py` - Complete Python simulation implementing Moran process models

### Research Output  
- `6999 Master's Project paper.pdf` - Full academic manuscript with methodology, results, and clinical implications

### Key Functions
- `TumorEvolutionModel` - Main simulation class with configurable parameters
- `moran_step()` - Single birth-death evolutionary event
- `calculate_diversity()` - Shannon entropy computation for genetic/phenotypic diversity
- `run_simulation()` - Complete evolutionary trajectory modeling

## Usage

### Basic Simulation
```python
# Initialize model
model = TumorEvolutionModel(genome_size=12, generations=200, fitness_type='additive')

# Run simulation  
genotypic_diversity, phenotypic_diversity, fitness_distributions = model.run_simulation()

# Analyze results
from scipy.stats import pearsonr
correlation, p_value = pearsonr(genotypic_diversity, phenotypic_diversity)
```

### Model Comparison
```python
# Compare paradigms
oncogene_results = plot_oncogene_addiction_models()
convergent_results = plot_convergent_evolution_model()

# Statistical analysis
calculate_statistics(results_dict)
```

## Model Limitations & Future Directions

### Current Limitations
- **Simplified fitness landscape** - Real cancer involves epistatic interactions
- **Binary genotype representation** - Actual mutations include CNVs, rearrangements, epigenetic changes
- **Uniform mutation rates** - Cancer genomes show mutational hotspots and context-dependency
- **Well-mixed populations** - Ignores spatial tumor architecture and microenvironmental effects

### Proposed Extensions
1. **Spatial modeling** - Incorporate tumor architecture and local selection pressures
2. **Epistatic fitness landscapes** - Model gene interaction effects
3. **Clinical validation** - Parameter estimation from real tumor evolution data
4. **Epigenetic inheritance** - Include non-genetic adaptive mechanisms

## Clinical Translation Potential

### Immediate Applications
- **Treatment selection algorithms** incorporating evolutionary potential assessments
- **Biomarker development** focusing on diversity metrics rather than static genetic markers
- **Clinical trial design** comparing evolution-guided vs. traditional therapeutic approaches

### Long-term Vision
Transition from static, gene-focused treatment paradigms to dynamic, evolution-based frameworks that account for tumor adaptive capacity and target multiple evolutionary escape routes simultaneously.

## Statistical Validation
The model successfully differentiates paradigms through quantitative metrics:
- **Correlation analysis** between genetic and phenotypic diversity
- **Shannon entropy** tracking over evolutionary time
- **Fitness distribution analysis** revealing convergent phenotypic outcomes

## Broader Significance
Beyond cancer, this framework applies to any system involving evolution under selection pressure: antimicrobial resistance, vaccine escape variants, and autoimmune disease progression. The principle that "driving species to extinction requires understanding all evolutionary pathways" has applications across medical and conservation biology.

## Academic Context
**Master's Project (BINF*6999)**  
University of Guelph, August 2025  
**Advisors**: Dr. Geoffrey Wood (Pathobiology), Dr. Ryan Gregory (Integrative Biology), Dr. Arijit Chakravarty (Fractal Therapeutics)
