module.exports = {
  preset: "ts-jest",
  collectCoverage: false,
  testEnvironment: "node",
  globals: {
    "ts-jest": {
      tsconfig: "tests/tsconfig.json",
      diagnostics: false,
    },
  },
  setupFilesAfterEnv: ["jest-expect-message"],
  transformIgnorePatterns: ["node_modules/(?!(ts3dutils)/)"],
  moduleNameMapper: {
    //"../src": "<rootDir>/lib/index.umd.min.js",
  },
}
