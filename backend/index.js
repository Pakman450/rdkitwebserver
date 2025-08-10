// index.js
require('dotenv').config();

const app = require('./app'); // import the actual app


const port = process.env.PORT || 3000;


console.log(`Running in ${process.env.NODE_ENV} mode, connecting to database: ${process.env.DB_NAME}`);

app.listen(port, () => {
  console.log(`Server running on http://localhost:${port}`);
});
