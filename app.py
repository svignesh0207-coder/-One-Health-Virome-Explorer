import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import entropy
import re

# =========================
# Page configuration
# =========================
st.set_page_config(
    page_title="One Health Virome Explorer",
    layout="wide"
)

# =========================
# Helper functions
# =========================

def extract_viral_family(taxon_name):
    n = taxon_name.lower()
    match = re.search(r"\b([a-z]+viridae)\b", n)
    return match.group(1).capitalize() if match else "Unresolved"


def classify_taxon_safe(taxon_name):
    n = taxon_name.lower()

    virus_type = (
        "Phage" if "phage" in n else
        "Eukaryotic virus" if "virus" in n else
        "Unknown"
    )

    family = extract_viral_family(taxon_name)

    if virus_type == "Phage":
        host, conf = "Bacterial", "High"
    elif any(x in n for x in ["herpes", "papilloma", "pox", "adeno"]):
        host, conf = "Mammal-associated", "High"
    elif any(x in n for x in ["avian", "gallid", "chicken", "fowl"]):
        host, conf = "Bird-associated", "Medium"
    elif any(x in n for x in ["baculovirus", "ascovirus", "nudivirus"]):
        host, conf = "Insect-associated", "High"
    else:
        host, conf = "Unknown", "Low"

    if host == "Mammal-associated":
        oh = "High"
    elif host in ["Bird-associated", "Insect-associated"]:
        oh = "Moderate"
    elif host == "Bacterial":
        oh = "Low"
    else:
        oh = "Uncertain"

    spill = (
        "Likely" if host == "Mammal-associated" else
        "Possible" if host in ["Bird-associated", "Insect-associated"] else
        "Not evident"
    )

    return family, host, conf, oh, spill


def calculate_alpha_diversity(counts):
    total = counts.sum()
    if total == 0:
        return 0, 0
    p = counts / total
    return entropy(p, base=np.e), 1 - np.sum(p ** 2)


# =========================
# App title
# =========================
st.title("🦠 One Health Virome Explorer")
st.caption(
    "Research-grade exploratory virome analysis integrating ecology, host inference, and One Health interpretation."
)

# =========================
# File upload
# =========================
uploaded_file = st.file_uploader(
    "Upload Kraken-style virome count table (CSV with columns: Taxon, Count)",
    type=["csv"]
)

if uploaded_file is None:
    st.info("⬆️ Upload a CSV file to begin analysis.")
    st.stop()

# =========================
# Load & validate data
# =========================
df = pd.read_csv(uploaded_file)

required_cols = {"Taxon", "Count"}
if not required_cols.issubset(df.columns):
    st.error("CSV must contain columns: Taxon and Count")
    st.stop()

df["Count"] = pd.to_numeric(df["Count"], errors="coerce").fillna(0).astype(int)
df = df[df["Count"] > 0].reset_index(drop=True)
total_reads = df["Count"].sum()

# =========================
# Classification
# =========================
df[[
    "Family_Assigned",
    "Host_Inference",
    "Host_Confidence",
    "OneHealth_Relevance",
    "Spillover_Potential"
]] = df["Taxon"].apply(lambda x: pd.Series(classify_taxon_safe(x)))

# =========================
# Diversity
# =========================
shannon, simpson = calculate_alpha_diversity(df["Count"])

# =========================
# Global One Health summary (used in multiple tabs)
# =========================
summary_lines = [
    f"Total viral taxa detected: {len(df)}",
    f"High One Health relevance taxa: {(df['OneHealth_Relevance'] == 'High').sum()}",
    f"Taxa with potential spillover relevance: {df['Spillover_Potential'].isin(['Likely', 'Possible']).sum()}",
    f"Environmental or unknown host taxa: {df['Host_Inference'].isin(['Environmental', 'Unknown']).sum()}",
    "",
    "Top 5 most abundant viral taxa:"
]

top5 = df.sort_values("Count", ascending=False).head(5)
for _, row in top5.iterrows():
    summary_lines.append(f"- {row['Taxon']} ({row['Count']} reads)")

one_health_summary_text = "\n".join(summary_lines)

# =========================
# Tabs
# =========================
tab0, tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
    "🏠 Overview",
    "📊 Community Structure",
    "🧬 Taxonomy & Host",
    "🌍 One Health Patterns",
    "🧫 Families",
    "🐾 Host-Specific",
    "🚨 Spillover View",
    "📥 Downloads"
])

# =========================
# TAB 0 — Overview
# =========================
with tab0:
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Viral Taxa", len(df))
    col2.metric("Total Reads", f"{total_reads:,}")
    col3.metric("Shannon Diversity", f"{shannon:.3f}")
    col4.metric("Simpson Diversity", f"{simpson:.3f}")

    st.markdown("---")
    st.subheader("🧠 Automated One Health Summary")
    st.text(one_health_summary_text)

    st.caption(
        "⚠️ Research use only. Inference-based host association and risk interpretation."
    )

# =========================
# TAB 1 — Community Structure
# =========================
with tab1:
    col1, col2 = st.columns(2)

    with col1:
        top_taxa = df.sort_values("Count", ascending=False).head(15)
        fig, ax = plt.subplots()
        sns.barplot(data=top_taxa, x="Count", y="Taxon", ax=ax)
        st.pyplot(fig)
        plt.close(fig)

    with col2:
        sorted_counts = df["Count"].sort_values(ascending=False).values
        ranks = np.arange(1, len(sorted_counts) + 1)
        fig, ax = plt.subplots()
        ax.plot(ranks, sorted_counts, marker="o", markersize=3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Rank (log)")
        ax.set_ylabel("Abundance (log)")
        st.pyplot(fig)
        plt.close(fig)

# =========================
# TAB 2 — Taxonomy & Host
# =========================
with tab2:
    host_filter = st.multiselect(
        "Host Inference",
        options=sorted(df["Host_Inference"].unique()),
        default=sorted(df["Host_Inference"].unique())
    )

    risk_filter = st.multiselect(
        "One Health Relevance",
        options=sorted(df["OneHealth_Relevance"].unique()),
        default=sorted(df["OneHealth_Relevance"].unique())
    )

    st.dataframe(
        df[
            df["Host_Inference"].isin(host_filter) &
            df["OneHealth_Relevance"].isin(risk_filter)
        ][[
            "Taxon", "Count", "Family_Assigned",
            "Host_Inference", "Host_Confidence",
            "OneHealth_Relevance", "Spillover_Potential"
        ]],
        use_container_width=True
    )

# =========================
# TAB 3 — One Health Patterns
# =========================
with tab3:
    col1, col2 = st.columns(2)

    with col1:
        fig, ax = plt.subplots()
        df["Host_Inference"].value_counts().plot.pie(
            autopct="%1.1f%%", ax=ax
        )
        st.pyplot(fig)
        plt.close(fig)

    with col2:
        fig, ax = plt.subplots()
        sns.barplot(
            x=df["OneHealth_Relevance"].value_counts().index,
            y=df["OneHealth_Relevance"].value_counts().values,
            ax=ax
        )
        st.pyplot(fig)
        plt.close(fig)

# =========================
# TAB 4 — Families
# =========================
with tab4:
    family_counts = (
        df[df["Family_Assigned"] != "Unresolved"]
        .groupby("Family_Assigned")["Count"]
        .sum()
        .sort_values(ascending=False)
        .head(10)
    )

    fig, ax = plt.subplots()
    sns.barplot(x=family_counts.values, y=family_counts.index, ax=ax)
    st.pyplot(fig)
    plt.close(fig)

# =========================
# TAB 5 — Host-Specific
# =========================
with tab5:
    selected_host = st.selectbox(
        "Select Host Group",
        sorted(df["Host_Inference"].unique())
    )

    host_df = df[df["Host_Inference"] == selected_host]
    st.dataframe(host_df.sort_values("Count", ascending=False).head(10))

    sh, si = calculate_alpha_diversity(host_df["Count"])
    col1, col2 = st.columns(2)
    col1.metric("Shannon", f"{sh:.3f}")
    col2.metric("Simpson", f"{si:.3f}")

# =========================
# TAB 6 — Spillover
# =========================
with tab6:
    spill_df = df[df["Spillover_Potential"].isin(["Likely", "Possible"])]

    if spill_df.empty:
        st.info("No spillover-relevant taxa detected.")
    else:
        st.dataframe(spill_df.sort_values("Count", ascending=False))

# =========================
# TAB 7 — Downloads
# =========================
with tab7:
    st.download_button(
        "Download Annotated Virome Table",
        df.to_csv(index=False),
        file_name="virome_annotated_table.csv"
    )

    st.download_button(
        "Download One Health Summary",
        one_health_summary_text,
        file_name="one_health_summary.txt"
    )
