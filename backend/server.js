// backend/server.js
const express = require('express');


const getFormula = require ('./services/pythonService.js')
const app = express();




app.get('/formula', async (req, res) => {
    const smiles = req.query.smiles || "CCO"; // ethanol

    try {
        const response = await getFormula.getFormula(smiles)
        res.json(response);
    } catch (err) {
        res.status(500).json({ error: err.message });
    }
});

app.listen(3000, () => console.log('Express listening on port 3000'));

module.exports = app
