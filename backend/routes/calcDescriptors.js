const express = require("express")
const router = express.Router()

const {
    calcDescriptors

}= require ('../services/pythonService.js')



router.get('/descriptors', async (req, res) => {

    const smiles = req.query.smiles || "CCCCCC"; // ethanol

    try {
        const response = await calcDescriptors(smiles)

        console.log(Object.keys(response))

        res.json(response);
    } catch (err) {
        res.status(500).json({ error: err.message });
    }
} 
)

module.exports = router