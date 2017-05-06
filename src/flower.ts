function makeFlower() {
	const golden = (1 + Math.sqrt(5)) / 2, rotateStep = 2 * PI * golden
	let flowerMesh = Mesh.plane().rotateZ(-45 * DEG).scale(1, 0.5, 1).rotateX(5 * DEG)
	const otherPetals = arrayFromFunction(200, i => {
		const f = 1 + 1/ 200 * i
		return flowerMesh
			.scale(f, f, f)
			.rotateY((i + 1) * 0.09 * DEG - 20 * DEG)
			.rotateZ(rotateStep * (i + 1))
	})
	flowerMesh = flowerMesh.concat.apply(flowerMesh, otherPetals)
	flowerMesh.compile()
	return flowerMesh
}
