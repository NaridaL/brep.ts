/**
 * Created by aval on 19.04.2017.
 */
abstract class FaceInfoFactory<T> {
	info(surface: Surface, contour: Edge[], holes: Edge[][]): T {
		throw new Error('no default implementation')
	}

	extrudeBottom(surface: Surface, contour: Edge[], holes: Edge[][] = []): T {
		return this.info(surface, contour, holes)
	}

	extrudeTop(surface: Surface, contour: Edge[], holes: Edge[][] = []): T {
		return this.info(surface, contour, holes)
	}

	extrudeWall(index: int, surface: Surface, contour: Edge[], holes: Edge[][] = []): T {
		return this.info(surface, contour, holes)
	}

	rotationWall(index: int, surface: Surface, contour: Edge[], holes: Edge[][] = []): T {
		return this.info(surface, contour, holes)
	}

	rotationStart(int, surface: Surface, contour: Edge[], holes: Edge[][] = []): T {
		return this.info(surface, contour, holes)
	}

	rotationEnd(surface: Surface, contour: Edge[], holes: Edge[][] = []): T {
		return this.info(surface, contour, holes)
	}

	newSubFace(original: Face, surface: Surface, contour: Edge[], holes: Edge[][] = []): T {
		return original.info
	}

	transform(original: Face, m4: M4, desc: string, surface: Surface, contour: Edge[], holes: Edge[][] = []): T {
		return original.info
	}
}