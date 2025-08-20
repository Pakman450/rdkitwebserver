// backend/server.js
const express = require('express');

const app = express();

const getFormula = require('./routes/getFormula')
const getAllDescriptors = require('./routes/calcDescriptors')
const getQED = require('./routes/getQED')

app.use('/v1', getFormula )

app.use('/v1', getAllDescriptors)

app.use('/v1', getQED)


module.exports = app
