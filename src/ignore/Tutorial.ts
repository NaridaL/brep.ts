/**
 * Created by aval on 08/02/2016.
 */
class Tutorial {
  readonly actions = []
  currentAction = -1

  constructor(readonly name: string, init) {
    this.name = name
    this.actions = []
    this.currentAction = -1
    init(this)
  }

  text(on, text) {
    this.actions.push({ type: "text", on: on, text: text })
  }

  wait(f) {
    this.actions.push({ type: "wait", f: f })
  }

  notifyEvent(e) {}

  startTutorial() {
    this.currentAction = 0
  }

  doText(action) {}
}

new Tutorial("Edge Champfer Tool", function (tut: Tutorial) {
  tut.text(
    undefined,
    `Save and click <a href="javascript:loadModel('edgeChampferTutorialModel')">here to load the required model`,
  )
  tut.wait(
    (e) => (e.type = "loadModel" && e.modelName == "edgeChampferTutorialModel"),
  )
  tut.text(undefined, "Select the highlighted edge")
  const edge1 = elByName("edge1")
  tut.setHighlight(edge1)
  tut.wait((e) => selected.equals([edge1]))
  tut.text("champferToolId", "Click on the Champfer Tool")
  tut.wait()
  tut.text(
    "champferEdgeSelect",
    "The selected edge shows up here, you can change it.",
  )
  tut.waitForContinue()
  tut.text("champferType", "...")
  tut.waitForContinue()
  tut.text(
    "champferOffset",
    "Adjust value to 3. Remember that you can use the scroll wheel.",
  )
})
