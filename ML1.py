import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import mean_absolute_error, accuracy_score, classification_report
import matplotlib.pyplot as plt

# Step 1: Parse the VCF File
def parse_vcf(file_path):
    data = []
    print("Reading VCF file...")

    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()

                if line.startswith("##") or line.startswith("#CHROM"):
                    continue  

                cols = line.split("\t")
                if len(cols) < 8:
                    print(f"Skipping malformed line: {line}")
                    continue

                chrom, pos, id_, ref, alt, qual, filter_, info, x, y, label = cols[:11]
                #label = cols[-1]
                try:
                    info_dict = {k: v for k, v in (item.split("=") for item in info.split(";") if "=" in item)}
                    dp = int(info_dict.get("DP", 0))
                    as_sb = float(info_dict.get("AS_SB_TABLE", 0).split("|")[0].split(',')[0])
                    pop_af = float(info_dict.get("POPAF", 0))

                    data.append({
                        "CHROM": chrom,
                        "POS": pos,
                        "ID": id_,
                        "REF": ref,
                        "ALT": alt,
                        "POP_AF": pop_af,
                        "DP": dp,
                        "AS_SB_TABLE": as_sb,
                    })
                except Exception as e:
                    print(f"Error parsing line: {e}")

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found!")
        return pd.DataFrame()  

    if not data:
        print("Error: No valid data parsed!")
        return pd.DataFrame()

    print(f"Successfully parsed {len(data)} variants!")
    return pd.DataFrame(data)

# Step 2: Normalize Data
def preprocess_data(df):
    if df.empty:
        print("Error: Cannot preprocess an empty DataFrame!")
        return df

    df = df.dropna()  

    for col in ["DP", "AS_SB_TABLE", "POP_AF"]:
        df[col] = (df[col] - df[col].min()) / (df[col].max() - df[col].min())

    print("Data normalized successfully!")
    return df

# Step 3: Label Variants
def label_variants(df, top_percent=10, bottom_percent=10):
    if df.empty:
        print("Error: Cannot label an empty DataFrame!")
        return df

    df = df.sort_values(by="POP_AF", ascending=False)  
    num_variants = len(df)

    top_cutoff = int(num_variants * top_percent / 100)
    bottom_cutoff = int(num_variants * bottom_percent / 100)

    df["Label"] = "Medium"
    df.iloc[:top_cutoff, df.columns.get_loc("Label")] = "High"  
    df.iloc[-bottom_cutoff:, df.columns.get_loc("Label")] = "Low"  

    print("Variant labels assigned successfully!")
    return df
def train_models(df):
    if df.empty:
        print("Error: Cannot train on an empty DataFrame!")
        return

    # Prepare features (X) and target labels (y)
    X = df[["DP", "AS_SB_TABLE", "POP_AF"]]  # Features
    y = df["Label"].map({"Low": 0, "Medium": 1, "High": 2})  # Convert labels to numbers

    # Split into train & test sets (80% training, 20% testing)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    print("Training data size:", X_train.shape)
    print("Testing data size:", X_test.shape)

    # Train Linear Regression model
    linear_model = LinearRegression()
    linear_model.fit(X_train, y_train)

    # Predict on the test set
    linear_preds = linear_model.predict(X_test)

    # Evaluate Linear Regression model
    mae = mean_absolute_error(y_test, linear_preds)
    print(f"Linear Regression Mean Absolute Error: {mae:.4f}")

    # Train Decision Tree model
    tree_model = DecisionTreeClassifier(random_state=42)
    tree_model.fit(X_train, y_train)

    # Predict on the test set
    tree_preds = tree_model.predict(X_test)

    # Evaluate Decision Tree model
    accuracy = accuracy_score(y_test, tree_preds)
    print(f"Decision Tree Accuracy: {accuracy:.4f}")

    # Print classification report
    print("Decision Tree Classification Report:")
    print(classification_report(y_test, tree_preds))

    # Plot Linear Regression Predictions
    plt.scatter(y_test, linear_preds, alpha=0.7)
    plt.xlabel("True Confidence Labels")
    plt.ylabel("Predicted Confidence Labels")
    plt.title("Linear Regression Predictions")
    plt.show()

    # Feature Importance in Decision Tree
    plt.barh(X.columns, tree_model.feature_importances_)
    plt.xlabel("Feature Importance")
    plt.title("Feature Importance (Decision Tree)")
    plt.show()

# Step 4: Main Function
def main():
    file_path = "/Users/neta/Documents/python_practice/mock_variants.vcf"  
    
    df = parse_vcf(file_path)

    if df.empty:
        print("No valid data was parsed. Exiting program.")
        return

    print("\nParsed VCF Data (First 5 rows):")
    print(df.head())

    df = preprocess_data(df)
    if df.empty:
        print("No valid data after preprocessing. Exiting program.")
        return

    print("\nNormalized Data (First 5 rows):")
    print(df.head())

    df = label_variants(df)
    if df.empty:
        print("No valid data after labeling. Exiting program.")
        return

    print("\nVariant Confidence Labels:")
    print(df["Label"].value_counts())

     # Train ML models
    train_models(df)

if __name__ == "__main__":
    main()

