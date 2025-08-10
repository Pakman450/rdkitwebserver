
//import 'dotenv/config'
require('dotenv').config();

async function getFormula(smiles) {
    const url = `${process.env.PYTHON_API_URL}/smiles_to_mol?smiles=${encodeURIComponent(smiles)}`;
    const res = await fetch(url);
    if (!res.ok) throw new Error(`Python API error: ${res.status}`);
    return res.json();
}

async function getQED(smiles) {
    const url = `${process.env.PYTHON_API_URL}/calc_qed?smiles=${encodeURIComponent(smiles)}`;
    const res = await fetch(url);
    if (!res.ok) throw new Error(`Python API error: ${res.status}`);
    return res.json();
}

module.exports = { getFormula , getQED }
