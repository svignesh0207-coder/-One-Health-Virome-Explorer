import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import entropy
import re
import io

# =========================
# Page configuration
# =========================
st.set_page_config(
    page_title="One Health Virome Explorer",
    layout="wide"
)

# =========================
# Custom CSS (Professional Theme)
# =========================
st.markdown("""
<style>
body {
    font-family: Arial, sans-serif;
}
h1, h2, h3 {
    color: #1f4e79;
}
</style>
""", unsafe_allow_html=True)

sns.set_style("whitegrid")

# =========================
# Helper functions
# =========================

def extract_viral_family(taxon):
    m = re.search(r"\b([a-z]+viridae)\b", taxon.lower())
    return m.group(1).capitalize() if m else "Unresolved"


def classify_taxon_safe(taxon):
    n = taxon.lower()

    if "phage" in n:
        host, conf = "Bacterial", "High"
    elif any(x in n for x in ["herpes", "papilloma", "pox", "adeno"]):
        host, conf = "Mammal-associated", "High"
    elif any(x in n for x in ["avian", "chicken", "fowl", "gallid"]):
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

    spill = "Likely" if host == "Mammal-associated" else \
            "Possible" if host in ["Bird-associated", "Insect-associated"] else \
            "Not evident"

    return extract_viral_family(taxon), host, conf, oh, spill


def calculate_diversity(counts):
    p = counts / counts.sum()
    return entropy(p), 1 - np.sum(p**2)


def fig_to_bytes(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    buf.seek(0)
    return buf

# =========================
# App Title
# =========================
st.title("🦠 One Health Virome Explorer")
st.caption(
    "Exploratory virome analysis integrating ecology, host inference, and One Health interpretation."
)

# =========================
# Upload Section
# =========================
uploaded = st.file_uploader(
    "Upload Kraken-style CSV (Taxon, Count)",
    type=["csv"]
)

if uploaded is None:
    st.info("⬆️ Upload a CSV file to begin analysis.")
    st.stop()

# =========================
# Data Loading (cached)
# =========================
with st.spinner("Processing virome data..."):
    df = pd.read_csv(uploaded)
    df["Count"] = pd.to_numeric(df["Count"], errors="coerce").fillna(0).astype(int)
    df = df[df["Count"] > 0].reset_index(drop=True)

    df[[
        "Family", "Host", "Host_Confidence",
        "OneHealth", "Spillover"
    ]] = df["Taxon"].apply(lambda x: pd.Series(classify_taxon_safe(x)))

total_reads = df["Count"].sum()
shannon, simpson = calculate_diversity(df["Count"])

# =========================
# Summary Text
# =========================
summary_text = "\n".join([
    f"Total viral taxa: {len(df)}",
    f"High One Health relevance: {(df['OneHealth']=='High').sum()}",
    f"Potential spillover taxa: {df['Spillover'].isin(['Likely','Possible']).sum()}",
    "",
    "Top 5 abundant taxa:",
    *[f"- {r.Taxon} ({r.Count})" for _, r in df.nlargest(5, "Count").iterrows()]
])

# =========================
# Tabs
# =========================
tabs = st.tabs([
    "🏠 Overview", "📊 Community", "🧬 Taxonomy",
    "🌍 One Health", "🧫 Families", "🐾 Hosts",
    "🚨 Spillover", "📥 Downloads"
])

# =========================
# TAB 0 — Overview
# =========================
with tabs[0]:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Taxa", len(df))
    c2.metric("Reads", f"{total_reads:,}")
    c3.metric("Shannon", f"{shannon:.3f}")
    c4.metric("Simpson", f"{simpson:.3f}")

    with st.expander("🧠 Automated One Health Summary"):
        st.text(summary_text)

# =========================
# TAB 1 — Community
# =========================
with tabs[1]:
    top = df.nlargest(15, "Count")

    fig, ax = plt.subplots()
    sns.barplot(data=top, x="Count", y="Taxon", ax=ax)
    st.pyplot(fig)
    st.download_button("⬇️ Download Plot", fig_to_bytes(fig), "top_taxa.png")
    plt.close(fig)

# =========================
# TAB 2 — Taxonomy
# =========================
with tabs[2]:
    search = st.text_input("🔍 Search Taxon")
    view = df[df["Taxon"].str.contains(search, case=False)] if search else df

    st.dataframe(
        view[[
            "Taxon", "Count", "Family",
            "Host", "Host_Confidence",
            "OneHealth", "Spillover"
        ]],
        use_container_width=True
    )

# =========================
# TAB 3 — One Health
# =========================
with tabs[3]:
    fig, ax = plt.subplots()
    df["OneHealth"].value_counts().plot.bar(ax=ax)
    st.pyplot(fig)
    plt.close(fig)

# =========================
# TAB 4 — Families
# =========================
with tabs[4]:
    fam = df.groupby("Family")["Count"].sum().nlargest(10)
    fig, ax = plt.subplots()
    sns.barplot(x=fam.values, y=fam.index, ax=ax)
    st.pyplot(fig)
    plt.close(fig)

# =========================
# TAB 5 — Hosts
# =========================
with tabs[5]:
    host = st.selectbox("Select Host", sorted(df["Host"].unique()))
    hdf = df[df["Host"] == host]
    st.dataframe(hdf.nlargest(10, "Count"))

# =========================
# TAB 6 — Spillover
# =========================
with tabs[6]:
    sdf = df[df["Spillover"].isin(["Likely", "Possible"])]
    if sdf.empty:
        st.info("No spillover-relevant taxa detected.")
    else:
        st.dataframe(sdf.nlargest(15, "Count"))

# =========================
# TAB 7 — Downloads
# =========================
with tabs[7]:
    st.download_button(
        "Download Annotated Table",
        df.to_csv(index=False),
        "virome_annotated.csv"
    )
    st.download_button(
        "Download One Health Summary",
        summary_text,
        "one_health_summary.txt"
    )
