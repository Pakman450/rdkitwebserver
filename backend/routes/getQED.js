const express = require("express")
const router = express.Router()

const {
    getQED

}= require ('../services/pythonService.js')


router.get('/getQED', async (req, res) => {
    const smiles = req.query.smiles || "CCO"; // ethanol

    try {
        const response = await getQED(smiles)

        console.log(res.body)

        res.json(response);
    } catch (err) {
        res.status(500).json({ error: err.message });
    }
} 
)

module.exports = router