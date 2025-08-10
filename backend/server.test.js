const request = require('supertest');
const app = require('./server');  // adjust the path if needed

describe('GET /formula', () => {
  it('returns formula for ethanol (default)', async () => {
    const res = await request(app).get('/formula');
    expect(res.statusCode).toBe(200);
    expect(res.body).toHaveProperty('formula');  // adjust if your response is different
  });

  it('returns formula for given SMILES', async () => {
    const res = await request(app).get('/formula').query({ smiles: 'CCO' });
    expect(res.statusCode).toBe(200);
    expect(res.body).toHaveProperty('formula');
    expect(res.body.formula).toBe('C2H6O');
  });

  
});
