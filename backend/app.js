// backend/server.js
const express = require('express');

const app = express();

const getFormula = require('./routes/getFormula')

app.use('/v1', getFormula )


module.exports = app
