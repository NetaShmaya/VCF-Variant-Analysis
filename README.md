# VCF Variant Analysis with Machine Learning

This project processes VCF (Variant Call Format) files to classify genetic variants into high-confidence, medium-confidence, and low-confidence categories using machine learning models. It uses Linear Regression and Decision Tree Classifier to predict variant confidence levels based on features like DP (Depth of Coverage), AS_SB_TABLE (Allele-Specific Strand Bias), and POP_AF (Population Allele Frequency).

Features
-Parses a VCF file and extracts relevant features
-Normalizes data for better model performance
-Labels variants into High, Medium, or Low confidence categories
-Trains machine learning models to classify variants
-Visualizes feature importance and model performance

Setup and Installation:

Clone this repository:
git clone <your-repo-url>
cd <your-repo-name>

(Optional) Create and activate a virtual environment:
python -m venv venv  
source venv/bin/activate  # On Windows: venv\Scripts\activate

Install dependencies:
pip install -r requirements.txt

Dependencies:
Python 3.x
pandas, scikit-learn, matplotlib