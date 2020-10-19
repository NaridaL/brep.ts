module.exports = {
  preset: "ts-jest",
  collectCoverage: false,
  testEnvironment: "node",
  globals: {
    "ts-jest": {
      tsConfig: {
        sourceMap: true,
        inlineSourceMap: true,
      },
      diagnostics: false,
    },
  },
  setupFilesAfterEnv: ["jest-expect-message"],
}
