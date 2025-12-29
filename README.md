🦠 One Health Virome Explorer

One Health Virome Explorer is an interactive Streamlit application for exploratory virome analysis of metagenomic outputs (e.g., Kraken2-style taxonomic count tables). The platform integrates ecological diversity metrics, rule-based host inference, viral family aggregation, and One Health–oriented interpretation to support hypothesis generation in virology, microbiome research, and environmental surveillance.

The tool is designed for research and exploratory use, emphasizing interpretability, transparency, and biological caution, rather than diagnostic or predictive claims.

⚠️ Research use only. This application is not intended for clinical diagnostics, regulatory decision-making, or direct epidemiological conclusions. All host associations, spillover relevance, and One Health interpretations are heuristic-based and should be validated using domain expertise and supporting evidence.

🎯 Key Highlights

Interdisciplinary Integration
Combines bioinformatics (taxonomic parsing), ecological analysis (Shannon & Simpson diversity), and One Health framing within a single interactive dashboard.

Scalable for Large Datasets
Efficiently handles large virome datasets (1,000+ taxa) using optimized Pandas operations and SciPy-based statistical computations.

Interpretation-Focused Design
Provides structured summaries, host-level views, and risk-oriented groupings that aid biological interpretation rather than raw visualization alone.

Extensible Architecture
Modular codebase designed to support future extensions such as bacterial modules, AMR integration, or model-based host prediction.

✨ Features
📤 Upload & Data Validation

Upload Kraken-style CSV files containing Taxon and Count columns

Automatic validation and filtering of non-zero counts

📊 Community Structure Analysis

Top taxa abundance plots

Rank–abundance curves (log–log scale) for dominance and evenness assessment

Alpha diversity metrics (Shannon entropy, Simpson index)

🧬 Taxonomy & Host Inference

Automated viral family extraction (regex-based from taxon names)

Rule-based host inference (bacterial, mammalian, avian, insect, unknown)

Confidence annotation for inferred host groups

🌍 One Health Patterns

Interactive filtering by host inference and One Health relevance

Distribution plots for host groups and relevance categories

🧫 Family-Level Insights

Aggregated abundance of dominant viral families

Rapid identification of family-level ecological patterns

🐾 Host-Specific Views

Drill-down into host-associated subsets

Host-specific dominant taxa and diversity metrics

🚨 Spillover Surveillance View

Dedicated view for taxa with Likely or Possible spillover relevance

Sorted by abundance to highlight dominant candidates

📥 Downloads & Reproducibility

Downloadable annotated virome tables (CSV)

Plain-text One Health summary for reporting and documentation

🛠️ Installation
Clone the Repository
git clone https://github.com/yourusername/one-health-virome-explorer.git
cd one-health-virome-explorer

(Recommended) Create a Virtual Environment
python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate

Install Dependencies
pip install -r requirements.txt


Required packages:

streamlit

pandas

numpy

matplotlib

seaborn

scipy

Run the App Locally
streamlit run app.py



📖 Usage

Upload Data
Provide a CSV file with at least:

Taxon (e.g., Escherichia phage RCS47)

Count (integer read or contig counts)

Navigate Tabs
Start with Overview for global metrics, then explore specialized tabs such as Community Structure, Families, and Spillover View.

Filter & Explore
Use interactive filters to subset taxa by host inference or One Health relevance.

Export Results
Download annotated tables and summaries for downstream analysis in R, Python, or manuscript preparation.

📄 Example Input CSV
Taxon,Count
Escherichia phage RCS47,4823
Salmonella phage S6Mu,4523
Hemileuca sp. nucleopolyhedrovirus,7

🧰 Technologies Used

Frontend: Streamlit

Data Processing: Pandas, NumPy

Visualization: Matplotlib, Seaborn

Statistics: SciPy (entropy-based diversity metrics)

Pattern Matching: Python regex for taxonomic parsing

This stack reflects best practices in applied bioinformatics tooling: modular, interpretable, and reproducible.

🤝 Contributing

Contributions from bioinformaticians, data scientists, and One Health researchers are welcome.

Suggestions include:

Multi-sample comparative analysis

Integration of bacterial and AMR modules

Model-based host prediction

Taxonomy enrichment via public databases

Please open an issue or submit a pull request.


🙏 Acknowledgments

Inspired by metagenomic visualization tools such as Kraken2 and Pavian, and developed to support exploratory research in environmental and poultry-associated viromics.

📬 Contact

Author: Vignesh S
GitHub:  https://github.com/svignesh0207-coder

LinkedIn:http://www.linkedin.com/in/vignesh-s-36b27b211

If you are a researcher or industry professional interested in One Health–oriented bioinformatics tools, feel free to connect.
