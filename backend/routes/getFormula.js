
const express = require("express")
const router = express.Router()

const {
    getFormula

}= require ('../services/pythonService.js')


router.get('/get_formula', async (req, res) => {
    const smiles = req.query.smiles || "CCO"; // ethanol

    try {
        const response = await getFormula(smiles)
        res.json(response);
    } catch (err) {
        res.status(500).json({ error: err.message });
    }
} 
)

module.exports = router