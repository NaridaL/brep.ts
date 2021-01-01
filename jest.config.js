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
  moduleNameMapper: {
    brepts: "<rootDir>/src/index.ts",
  },
}
