const request = require('supertest');
const app = require('./app');  // adjust the path if needed

describe('GET /v1/formula', () => {
  it('returns formula for ethanol (default)', async () => {
    const res = await request(app).get('/v1/get_formula');
    expect(res.statusCode).toBe(200);
    expect(res.body).toHaveProperty('formula');  // adjust if your response is different
  });

  it('returns formula for given SMILES', async () => {
    const res = await request(app).get('/v1/get_formula').query({ smiles: 'CCO' });
    expect(res.statusCode).toBe(200);
    expect(res.body).toHaveProperty('formula');
    expect(res.body.formula).toBe('C2H6O');
  });


});
