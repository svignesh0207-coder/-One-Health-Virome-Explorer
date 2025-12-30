import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from scipy.stats import entropy
import re
import io
import zipfile

# =========================
# Page config
# =========================
st.set_page_config(
    page_title="One Health Virome Explorer",
    layout="wide"
)

# =========================
# Helper functions
# =========================

def extract_family(taxon):
    m = re.search(r"\b([a-z]+viridae)\b", taxon.lower())
    return m.group(1).capitalize() if m else "Unresolved"


def classify_taxon(taxon):
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

    return extract_family(taxon), host, conf, oh, spill


def alpha_diversity(counts):
    p = counts / counts.sum()
    return entropy(p), 1 - np.sum(p**2)


# =========================
# Title
# =========================
st.title("🦠 One Health Virome Explorer")
st.caption("Interactive exploratory virome analysis with One Health framing.")

# =========================
# Upload
# =========================
uploaded = st.file_uploader(
    "Upload Kraken-style CSV (Taxon, Count)",
    type=["csv"]
)

if uploaded is None:
    st.info("⬆️ Upload a CSV file to begin.")
    st.stop()

# =========================
# Load data
# =========================
with st.spinner("Processing virome data..."):
    df = pd.read_csv(uploaded)
    df["Count"] = pd.to_numeric(df["Count"], errors="coerce").fillna(0).astype(int)
    df = df[df["Count"] > 0].reset_index(drop=True)

    df[["Family", "Host", "Host_Confidence", "OneHealth", "Spillover"]] = (
        df["Taxon"].apply(lambda x: pd.Series(classify_taxon(x)))
    )

total_reads = df["Count"].sum()
shannon, simpson = alpha_diversity(df["Count"])

# =========================
# Summary
# =========================
summary_text = "\n".join([
    f"Total viral taxa: {len(df)}",
    f"High One Health relevance taxa: {(df['OneHealth']=='High').sum()}",
    f"Potential spillover taxa: {df['Spillover'].isin(['Likely','Possible']).sum()}",
    "",
    "Top 5 abundant taxa:",
    *[f"- {r.Taxon} ({r.Count})" for _, r in df.nlargest(5, 'Count').iterrows()]
])

# =========================
# Tabs
# =========================
tabs = st.tabs([
    "🏠 Overview", "📊 Community",
    "🧬 Taxonomy", "🌍 One Health",
    "🧫 Families", "🐾 Hosts",
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
# TAB 1 — Community Structure
# =========================
with tabs[1]:
    top_n = st.slider("Top N taxa", 5, 30, 15)
    top_df = df.nlargest(top_n, "Count")

    fig_bar = px.bar(
        top_df,
        x="Count",
        y="Taxon",
        orientation="h",
        title="Top Viral Taxa by Abundance"
    )
    st.plotly_chart(fig_bar, use_container_width=True)

    sorted_counts = df["Count"].sort_values(ascending=False).values
    ranks = np.arange(1, len(sorted_counts) + 1)

    fig_rank = px.line(
        x=ranks,
        y=sorted_counts,
        log_x=True,
        log_y=True,
        labels={"x": "Rank", "y": "Abundance"},
        title="Rank–Abundance Curve"
    )
    st.plotly_chart(fig_rank, use_container_width=True)

# =========================
# TAB 2 — Taxonomy
# =========================
with tabs[2]:
    search = st.text_input("🔍 Search taxon")
    view = df[df["Taxon"].str.contains(search, case=False)] if search else df

    st.dataframe(
        view[[
            "Taxon", "Count", "Family",
            "Host", "Host_Confidence",
            "OneHealth", "Spillover"
        ]],
        use_container_width=True,
        hide_index=True
    )

# =========================
# TAB 3 — One Health Patterns
# =========================
with tabs[3]:
    col1, col2 = st.columns(2)

    # --- Host composition pie chart ---
    with col1:
        fig_pie = px.pie(
            df,
            names="Host",
            title="Host Group Composition"
        )
        st.plotly_chart(fig_pie, use_container_width=True)

    # --- One Health relevance bar chart (SAFE FIX) ---
    with col2:
        oh_counts = (
            df["OneHealth"]
            .value_counts()
            .reset_index()
        )
        oh_counts.columns = ["OneHealth", "Count"]  # 🔑 CRITICAL LINE

        fig_oh = px.bar(
            oh_counts,
            x="OneHealth",
            y="Count",
            title="One Health Relevance Distribution",
            labels={
                "OneHealth": "One Health Category",
                "Count": "Number of Viral Taxa"
            }
        )
        st.plotly_chart(fig_oh, use_container_width=True)
# =========================
# TAB 4 — Families
# =========================
with tabs[4]:
    fam_df = (
        df[df["Family"] != "Unresolved"]
        .groupby("Family")["Count"]
        .sum()
        .reset_index()
        .nlargest(10, "Count")
    )

    fig_fam = px.bar(
        fam_df,
        x="Count",
        y="Family",
        orientation="h",
        title="Top Viral Families"
    )
    st.plotly_chart(fig_fam, use_container_width=True)

# =========================
# TAB 5 — Host-specific
# =========================
with tabs[5]:
    host = st.selectbox("Select Host Group", sorted(df["Host"].unique()))
    hdf = df[df["Host"] == host]

    st.dataframe(
        hdf.nlargest(10, "Count"),
        use_container_width=True,
        hide_index=True
    )

# =========================
# TAB 6 — Spillover
# =========================
with tabs[6]:
    sdf = df[df["Spillover"].isin(["Likely", "Possible"])]

    if sdf.empty:
        st.info("No spillover-relevant taxa detected.")
    else:
        st.dataframe(
            sdf.nlargest(15, "Count"),
            use_container_width=True,
            hide_index=True
        )

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

    # ZIP all visuals
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w") as z:
        z.writestr("top_taxa.html", fig_bar.to_html())
        z.writestr("rank_abundance.html", fig_rank.to_html())
        z.writestr("host_pie.html", fig_pie.to_html())
        z.writestr("one_health_bar.html", fig_oh.to_html())
        z.writestr("family_bar.html", fig_fam.to_html())

    zip_buffer.seek(0)

    st.download_button(
        "⬇️ Download ALL Visuals (ZIP)",
        zip_buffer,
        "one_health_virome_visuals.zip",
        mime="application/zip"
    )



