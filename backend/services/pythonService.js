
import 'dotenv/config'


export async function getFormula(smiles) {
    const url = `${process.env.PYTHON_API_URL}/smiles-to-mol?smiles=${encodeURIComponent(smiles)}`;
    const res = await fetch(url);
    if (!res.ok) throw new Error(`Python API error: ${res.status}`);
    return res.json();
}

