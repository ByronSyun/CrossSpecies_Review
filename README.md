# Cross-Species Transfer Learning for Splicing Variant Prediction

This repository contains the code for evaluating splicing variant prediction models across multiple species, as described in our manuscript:

**"Cross-Species Transfer Learning for Splicing Variant Prediction: A Review of Task-Specific Deep Learning and Genomic Foundation Models"**

## 📋 Repository Structure

The repository is organized into **4 main directories**, one for each species/dataset:

```
.
├── splicevardb/          # Human benchmark (SpliceVarDB)
├── ratGTEx/              # Rat benchmark (RatGTEx)
├── pigGTEx/              # Pig benchmark (PigGTEx)
└── chickenGTEx/          # Chicken benchmark (ChickenGTEx)
```

Each species directory follows the same structure:

```
<species>/
├── preprocessing/        # Data preprocessing scripts
│   ├── 01_create_positive_samples.py
│   ├── 02_create_negative_pool.py
│   ├── 03_create_balanced_benchmark.py
│   └── ...
├── baselines/           # Model evaluation scripts
│   ├── AlphaGenome/
│   ├── Evo2/
│   ├── DNABERT-2/
│   ├── SpliceAI/
│   ├── Pangolin/
│   ├── MMSplice/
│   ├── SpliceTransformer/
│   ├── SpliceBERT/
│   ├── Nucleotide_Transformer/
│   ├── plot_<species>_performance.py    # ROC/PRC curves
│   └── visualisation/                   # Output figures
└── README.md (optional)
```

---

## 🚀 Quick Start

### 1. Environment Setup

Each model has its own environment requirements. Please refer to the `README.md` in each model's directory (e.g., `ratGTEx/baselines/Evo2/README.md`).

### 2. Data Preprocessing

```bash
# Example: Rat benchmark preprocessing
cd ratGTEx/preprocessing
python 01_create_positive_samples.py
python 02_create_negative_pool.py
python 03_create_balanced_benchmark.py
```

### 3. Run Model Evaluations

```bash
# Example: Run Evo2 on Rat
cd ratGTEx/baselines/Evo2
bash run_evo2_rat_server.sh
```

### 4. Generate Figures

```bash
# Example: Generate ROC/PRC curves for Rat
cd ratGTEx/baselines
python plot_ratgtex_performance.py

# Generate cross-species comparison (Figure 6)
python plot_cross_species_comparison.py
```

---

## 📊 Models Evaluated

### Task-Specific Deep Learning Models
- **SpliceAI**: CNN-based splice site predictor
- **Pangolin**: Dilated CNN with evolutionary features
- **MMSplice**: Modular, interpretable framework
- **SpliceTransformer**: Context-aware Transformer

### Genomic Foundation Models (GFMs)
- **AlphaGenome**: Supervised multi-task GFM
- **Evo2**: 7B-parameter generative GFM (StripedHyena architecture)
- **DNABERT-2**: 117M-parameter masked language model
- **Nucleotide Transformer**: 500M-2.5B parameter Transformer
- **SpliceBERT**: Splice-specialized pre-trained model

---

## 📁 Data Files

**Note**: This repository contains **only code** (Python scripts, shell scripts, R scripts, and documentation). **Data files are excluded** due to their large size and are governed by the following `.gitignore` rules:

```gitignore
# Excluded data files
*.csv
*.tsv
*.json
*.jsonl
*.npz
*.npy
*.h5
*.hdf5
*.vcf
*.vcf.gz
*.pkl
*.pickle
*.joblib
*.fasta
*.fa
```

### Where to Get Data

- **Human (SpliceVarDB)**: Download from [https://compbio.ccia.org.au/splicevardb/](https://compbio.ccia.org.au/splicevardb/)
- **Rat (RatGTEx)**: Contact the authors or download from [RatGTEx portal](https://ratgtex.org/)
- **Pig (PigGTEx)**: Part of FarmGTEx consortium data
- **Chicken (ChickenGTEx)**: Part of FarmGTEx consortium data

---

## 📖 Citation

If you use this code in your research, please cite our paper:

```bibtex
@article{CrossSpeciesReview2025,
  title={Cross-Species Transfer Learning for Splicing Variant Prediction: A Review of Task-Specific Deep Learning and Genomic Foundation Models},
  author={[Authors]},
  journal={Briefings in Bioinformatics},
  year={2025},
  note={In preparation}
}
```

---

## 🔧 Key Scripts

### Preprocessing
- `<species>/preprocessing/01_create_positive_samples.py`: Extract sQTL-derived positive samples
- `<species>/preprocessing/02_create_negative_pool.py`: Generate high-quality negative controls
- `<species>/preprocessing/03_create_balanced_benchmark.py`: Create balanced benchmark datasets

### Model Evaluation
- `<species>/baselines/<model>/predict_<species>_<model>.py`: Run model predictions
- `<species>/baselines/<model>/evaluate_<model>_<species>.py`: Evaluate model performance
- `<species>/baselines/plot_<species>_performance.py`: Generate ROC/PRC curves

### Cross-Species Comparison
- `ratGTEx/baselines/plot_cross_species_comparison.py`: Generate Figure 6 (cross-species bar chart)

---

## 📝 Important Notes

### Figure 5 (ROC/PRC Curves)
- **MMSplice is excluded** from cross-species visualization due to catastrophic coverage collapse (5-12%)
- Models shown: AlphaGenome, Evo2, Pangolin, SpliceAI, SpliceTransformer

### Figure 6 (Cross-Species Comparison Bar Chart)
- **MMSplice is retained** to demonstrate Human → Cross-species performance degradation
- Models shown: All 6 models including MMSplice

---

## 🛠️ Requirements

- Python 3.8+
- PyTorch 2.0+
- Transformers 4.30+
- scikit-learn 1.2+
- pandas, numpy, matplotlib, seaborn
- R 4.0+ (for Sankey diagrams)

Specific model requirements are detailed in each model's subdirectory.

---

## 📧 Contact

For questions about the code or data, please contact:
- **Byron Sun**: [your.email@example.com]

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- SpliceVarDB team for providing the human benchmark dataset
- FarmGTEx consortium for pig and chicken data
- RatGTEx team for rat data
- All model authors for making their code publicly available

