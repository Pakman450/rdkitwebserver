
//import 'dotenv/config'
require('dotenv').config();

async function getFormula(smiles) {
    const url = `${process.env.PYTHON_API_URL}/smiles_to_mol?smiles=${encodeURIComponent(smiles)}`;
    const res = await fetch(url);
    if (!res.ok) throw new Error(`Python API error: ${res.status}`);
    return res.json();
}

async function getQED(smiles) {
    const url = `${process.env.PYTHON_API_URL}/v1/calc_qed?smiles=${encodeURIComponent(smiles)}`;
    const res = await fetch(url);
    if (!res.ok) throw new Error(`Python API error: ${res.status}`);
    return res.json();
}

async function calcDescriptors(smiles) {

    const query = smiles
        .map(s => `smiles=${encodeURIComponent(s)}`)
        .join('&')

    console.log(query)

    const url = `${process.env.PYTHON_API_URL}/v1/descriptors?${query}`;
    const res = await fetch(url);
    if (!res.ok) throw new Error(`Python API error: ${res.status}`);
    return res.json();
}


module.exports = { getFormula , getQED, calcDescriptors }
