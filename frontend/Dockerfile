FROM node:18-alpine

# Set working directory
WORKDIR /app

# Install dependencies
COPY package.json package-lock.json ./
RUN npm install
RUN npm install react-scripts
RUN npm install --save-dev @babel/plugin-proposal-private-property-in-object

# Copy project files
COPY . .

# Debug: list files and contents
RUN ls -la
RUN ls -la src
RUN ls -la public

# Build the app with verbose output
RUN npm run build --verbose

# Debug: list build contents
RUN ls -la build

# Install serve
RUN npm install -g serve

# Expose port
EXPOSE 3000

# Start the app in development mode
CMD ["npm", "start"]
