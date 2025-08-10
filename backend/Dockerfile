# Start from official Node image
FROM node:20

# Set working directory
WORKDIR /app

# Copy package files and install deps
COPY package*.json ./
RUN npm install --production

# Copy backend code
COPY . .

# Expose Express port
EXPOSE 3000

# Start Express server
CMD ["node", "server.js"]
