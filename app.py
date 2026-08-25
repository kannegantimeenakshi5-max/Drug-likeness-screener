import pandas as pd
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

st.set_page_config(page_title="Cheminformatics Drug Screener", layout="wide")
st.title("💊 Cheminformatics Lead Screener & Lipinski Analyzer")

def calculate_lipinski(smiles_list):
    results = []
    for sm in smiles_list:
        mol = Chem.MolFromSmiles(sm)
        if mol is None:
            continue
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        
        # Lipinski Violations
        violations = sum([mw > 500, logp > 5.0, hbd > 5, hba > 10])
        passed = violations <= 1
        
        results.append({
            "SMILES": sm,
            "MW (Da)": round(mw, 2),
            "LogP": round(logp, 2),
            "HBD": hbd,
            "HBA": hba,
            "Rotatable_Bonds": rot_bonds,
            "Violations": violations,
            "Lipinski_Pass": passed
        })
    return pd.DataFrame(results)

# Sample SMILES data
sample_smiles = "CC(=O)Oc1ccccc1C(=O)O\nCC(=O)Nc1ccc(O)cc1\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C\nCC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"

user_input = st.sidebar.text_area("Enter SMILES strings (one per line):", value=sample_smiles, height=150)
smiles_list = [s.strip() for s in user_input.split("\n") if s.strip()]

if st.sidebar.button("Run Screening Pipeline"):
    df = calculate_lipinski(smiles_list)
    
    st.subheader("📊 Molecular Properties & Filter Results")
    st.dataframe(df, use_container_width=True)
    
    col1, col2 = st.columns(2)
    with col1:
        fig, ax = plt.subplots()
        sns.scatterplot(data=df, x="MW (Da)", y="LogP", hue="Lipinski_Pass", palette="coolwarm", s=100, ax=ax)
        ax.axhline(5, color='grey', linestyle='--', alpha=0.7)
        ax.axvline(500, color='grey', linestyle='--', alpha=0.7)
        ax.set_title("LogP vs Molecular Weight")
        st.pyplot(fig)
        
    with col2:
        fig2, ax2 = plt.subplots()
        sns.barplot(data=df, x=df.index, y="Violations", palette="viridis", ax=ax2)
        ax2.set_title("Rule Violations per Molecule")
        ax2.set_xlabel("Compound Index")
        st.pyplot(fig2)
        
    csv = df.to_csv(index=False).encode('utf-8')
    st.download_button("Download Filtered Leads (CSV)", data=csv, file_name="screened_leads.csv", mime="text/csv")
