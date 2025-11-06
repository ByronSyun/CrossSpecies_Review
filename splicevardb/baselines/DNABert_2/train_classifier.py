import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score, confusion_matrix
import joblib
import argparse
import seaborn as sns
import matplotlib.pyplot as plt

def train_and_evaluate(args):
    """
    Loads data, trains a logistic regression classifier, evaluates its performance,
    and saves the trained model.
    """
    # Load the pre-generated embeddings and labels
    print(f"Loading data from {args.input_file}...")
    data = np.load(args.input_file, allow_pickle=True)
    embeddings = data['embeddings']
    labels = data['labels']
    variant_ids = data['variant_ids']
    
    print(f"Loaded {len(embeddings)} samples.")
    
    # Split the data into training and testing sets (80/20 split)
    # Using stratify to ensure the proportion of labels is the same in train and test sets
    X_train, X_test, y_train, y_test, _, variant_ids_test = train_test_split(
        embeddings, labels, variant_ids, test_size=0.2, random_state=42, stratify=labels
    )
    
    print(f"Training set size: {len(X_train)}")
    print(f"Test set size: {len(X_test)}")
    
    # Initialize and train the Logistic Regression model
    print("Training Logistic Regression model...")
    model = LogisticRegression(random_state=42, max_iter=1000, class_weight='balanced')
    model.fit(X_train, y_train)
    
    # Make predictions on the test set
    print("Evaluating model performance on the test set...")
    y_pred = model.predict(X_test)
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    
    # Evaluate the model
    accuracy = accuracy_score(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_pred_proba)
    
    print("\n--- Model Evaluation Results ---")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"AUC-ROC Score: {roc_auc:.4f}")
    
    print("\nClassification Report:")
    # target_names=['normal', 'splice-altering']
    print(classification_report(y_test, y_pred, target_names=['normal (0)', 'splice-altering (1)']))
    
    print("\nConfusion Matrix:")
    cm = confusion_matrix(y_test, y_pred)
    print(cm)
    
    # Plotting the confusion matrix
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=['Predicted Normal', 'Predicted Splice-Altering'],
                yticklabels=['Actual Normal', 'Actual Splice-Altering'])
    plt.title('Confusion Matrix')
    plt.ylabel('Actual Label')
    plt.xlabel('Predicted Label')
    plt.savefig(args.confusion_matrix_file)
    print(f"Confusion matrix plot saved to {args.confusion_matrix_file}")
    
    # Save the trained model to disk
    print(f"\nSaving trained model to {args.model_output_file}...")
    joblib.dump(model, args.model_output_file)
    print("Model saved successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train and evaluate a classifier on DNABERT-2 embeddings.")
    parser.add_argument(
        "--input_file", 
        type=str, 
        default="dnabert2_predictions.npz",
        help="Path to the input NPZ file containing embeddings and labels."
    )
    parser.add_argument(
        "--model_output_file", 
        type=str, 
        default="dnabert2_classifier.joblib",
        help="Path to save the trained classifier model."
    )
    parser.add_argument(
        "--confusion_matrix_file", 
        type=str, 
        default="confusion_matrix.png",
        help="Path to save the confusion matrix plot."
    )
    
    args = parser.parse_args()
    train_and_evaluate(args) 