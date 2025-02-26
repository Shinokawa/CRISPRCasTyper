# CRISPR-Cas Predictor

This project is designed to analyze protein sequences and predict potential CRISPR-Cas proteins using various computational methods. It integrates multiple modules for sequence analysis, HMMER searches, and classification of CRISPR-Cas proteins.

## Project Structure

```
crispr-cas-predictor
├── src
│   ├── __init__.py
│   ├── main.py               # Entry point for the application
│   ├── protein_analyzer.py   # Contains ProteinAnalyzer class for sequence analysis
│   ├── hmmer_search.py        # Contains HMMERSearch class for HMMER tool integration
│   ├── cas_typing.py          # Contains CASTyping class for CRISPR-Cas protein classification
│   ├── utils
│   │   ├── __init__.py
│   │   ├── file_handling.py   # Utility functions for file handling
│   │   └── sequence_tools.py   # Utility functions for sequence processing
│   └── models
│       ├── __init__.py
│       └── cas_classifier.py   # Contains CASClassifier class for machine learning model
├── data
│   ├── hmm_models
│   │   └── cas_profiles.hmm    # HMM model file for HMMER searches
│   └── example_sequences
│       └── test_proteins.fasta  # Example protein sequences for testing
├── tests
│   ├── __init__.py
│   ├── test_hmmer_search.py     # Unit tests for HMMER search functionality
│   └── test_cas_typing.py       # Unit tests for CRISPR-Cas typing functionality
├── notebooks
│   └── analysis_examples.ipynb   # Jupyter Notebook for analysis examples and visualizations
├── requirements.txt              # List of required Python packages
├── setup.py                      # Installation script with package metadata
└── README.md                     # Project documentation and usage instructions
```

## Installation

To install the required dependencies, run:

```
pip install -r requirements.txt
```

## Usage

To run the application, use the following command:

```
python src/main.py <input_fasta_file> <output_directory>
```

Replace `<input_fasta_file>` with the path to your input FASTA file containing protein sequences, and `<output_directory>` with the desired output directory for results.

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue for any suggestions or improvements.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.