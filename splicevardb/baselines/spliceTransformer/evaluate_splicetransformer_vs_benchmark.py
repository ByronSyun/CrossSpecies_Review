import argparse
from typing import Optional
import pandas as pd
from sklearn.metrics import (
    roc_auc_score,
    precision_recall_curve,
    auc,
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
)


def _normalise_chrom(chrom: str) -> str:
    chrom = str(chrom)
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def _build_variant_id(chrom: str, pos_1based: int, ref: str, alt: str) -> str:
    return f"{_normalise_chrom(chrom)}-{int(pos_1based)}-{ref}-{alt}"


def read_predictions(pred_path: str) -> pd.DataFrame:
    # SpliceTransformer output is CSV-like with leading index column
    try:
        df = pd.read_csv(pred_path, sep=None, engine="python")
    except Exception:
        df = pd.read_csv(pred_path)

    # Handle possible unnamed index column
    if df.columns[0].startswith("Unnamed") or df.columns[0] == "":
        df = df.drop(columns=[df.columns[0]])

    # Minimal required columns
    expected = {"#CHROM", "POS", "REF", "ALT", "score"}
    missing = expected.difference(df.columns)
    if missing:
        raise ValueError(f"Predictions file missing columns: {missing}")

    df["variant_id"] = [
        _build_variant_id(c, p, r, a)
        for c, p, r, a in zip(df["#CHROM"], df["POS"], df["REF"], df["ALT"])
    ]
    df = df[["variant_id", "score"]]
    return df


def read_benchmark(bench_path: str) -> pd.DataFrame:
    df = pd.read_csv(bench_path, sep="\t", dtype=str)
    # Expect columns: hg38, classification, location (at least)
    needed = {"hg38", "classification"}
    if not needed.issubset(df.columns):
        raise ValueError(f"Benchmark file must contain columns {needed}, got {list(df.columns)}")

    def parse_hg38(x: str):
        parts = str(x).strip('"').split("-")
        if len(parts) != 4:
            return None
        chrom, pos, ref, alt = parts
        chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
        try:
            pos_1 = int(pos)
        except Exception:
            return None
        return chrom, pos_1, ref, alt

    parsed = df["hg38"].apply(parse_hg38)
    mask = parsed.notnull()
    df = df[mask].copy()
    parsed = parsed[mask]
    df[["chrom", "pos_1based", "ref", "alt"]] = pd.DataFrame(parsed.tolist(), index=df.index)

    df["label"] = df["classification"].str.strip().str.lower().map(
        {"splice-altering": 1, "splicealtering": 1, "normal": 0}
    )
    df = df.dropna(subset=["label"]).copy()
    df["label"] = df["label"].astype(int)

    df["variant_id"] = [
        _build_variant_id(c, p, r, a)
        for c, p, r, a in zip(df["chrom"], df["pos_1based"], df["ref"], df["alt"])
    ]
    return df[["variant_id", "label"]]


def evaluate(pred_path: str, bench_path: str, out_path: Optional[str] = None, threshold: float = 0.2) -> None:
    pred_df = read_predictions(pred_path)
    bench_df = read_benchmark(bench_path)

    merged = pd.merge(pred_df, bench_df, on="variant_id", how="inner")
    if merged.empty:
        raise RuntimeError("No overlapping variants between predictions and benchmark.")

    y_true = merged["label"]
    y_scores = merged["score"]
    y_pred = (y_scores >= threshold).astype(int)

    auroc = roc_auc_score(y_true, y_scores)
    pr, rc, _ = precision_recall_curve(y_true, y_scores)
    auprc = auc(rc, pr)
    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred)
    rec = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)
    cm = confusion_matrix(y_true, y_pred)

    print(f"Merged: {len(merged)}")
    print(f"AUROC: {auroc:.4f} | AUPRC: {auprc:.4f}")
    print(f"Accuracy: {acc:.4f} | Precision: {prec:.4f} | Recall: {rec:.4f} | F1: {f1:.4f}")
    print("Confusion matrix:\n", cm)

    if out_path:
        merged.to_csv(out_path, sep="\t", index=False)
        print(f"Saved merged results -> {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate SpliceTransformer predictions vs SpliceVarDB benchmark.")
    parser.add_argument("--pred", required=True, help="Path to SpliceTransformer predictions file (CSV-like with score column).")
    parser.add_argument("--benchmark", required=True, help="Path to benchmark TSV (splicevardb_filter_benchmark.tsv).")
    parser.add_argument("--out", help="Optional path to save merged evaluation table.")
    parser.add_argument("--threshold", type=float, default=0.2, help="Decision threshold for classification metrics.")
    args = parser.parse_args()

    evaluate(args.pred, args.benchmark, args.out, args.threshold)


