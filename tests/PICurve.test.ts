import { suite, test, testCurve, testISTs, testPointT } from './manager'

import { DEG, V, V3, M4, eq } from 'ts3dutils'
import {
	BezierCurve,
	EllipsoidSurface,
	P3,
	PICurve,
	PlaneSurface,
	ProjectedCurveSurface,
	SemiEllipsoidSurface,
	SemiCylinderSurface,
	SemiEllipseCurve,
} from '..'

suite('PICurve', () => {
	suite('pointT', () => {
		test('1', assert => {
			const pic = PICurve.forStartEnd(
				new ProjectedCurveSurface(
					new BezierCurve(
						V(0.20296874999999998, 0.11703124999999998, 1.2),
						V(0.20296874999999998, 0.2240625, 1.2),
						V(0.14500000000000002, 0.2890625, 1.2),
						V(0.010937499999999989, 0.2890625, 1.2),
						0,
						1,
					),
					V(0, 0, -1),
					0,
					1,
					-Infinity,
					Infinity,
				),
				new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z),
				V(0.010937499999999989, 0.2890625, -0.9572477433702835),
				V(0.20296874999999998, 0.11703124999999998, -0.9721663299286162),
				0.02,
			)
			const p = V(0.010937499999999989, 0.2890625, 0.9572477433702835)
			testPointT(assert, pic, p, NaN)
		})
		test('2', assert => {
			const pic = PICurve.forParametricStartEnd(
				new ProjectedCurveSurface(
					new BezierCurve(
						V(0, 2, 0),
						V(0.3333333333333333, 1, 0),
						V(0.6666666666666666, -1, 0),
						V(1, -2, 0),
						-0.1,
						1.1,
					),
					V(0, 0, -1),
					0,
					1,
					-1,
					0,
				),
				new SemiCylinderSurface(
					new SemiEllipseCurve(
						V(0, 0, 0.5),
						V(1.2246467991473533e-17, 0, -0.2),
						V(0, 0.2, 0),
						0,
						3.141592653589793,
					),
					V(1, 0, 6.123233995736766e-17),
					0,
					3.141592653589793,
					0,
					2,
				),
				V(0.543312790740714, -0.5455452074045042, 0),
				V(0.45794084603336055, -0.4349260970573627, 0),
				0.05,
				V(0.0026015534861132676, 0.04993227332556462, 0),
				0,
				10,
			)
			const p = pic.points[8]
			testPointT(assert, pic, p, 8)
		})
		test('3', assert => {
			const pic = PICurve.forParametricStartEnd(
				new ProjectedCurveSurface(
					new BezierCurve(
						V(0.01100000000000001, 0.28900000000000003, 1.2),
						V(-0.04100000000000001, 0.28900000000000003, 1.2),
						V(-0.1, 0.279, 1.2),
						V(-0.16499999999999998, 0.255, 1.2),
						0,
						1,
					),
					V(0, 0, -1),
					0,
					1,
					0,
					1,
				),
				new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z),
				V(-0.19976905595802608, 0.24321160282169424, 0),
				V(1.0002090044741139, 0.24724388010624657, 0),
				0.05,
				V(0.049999788672478875, -0.000145370930566841, 0),
				0,
				24,
			)
			const p = V(-0.16499999999999998, 0.255, 0.9527591510974849)
			testPointT(assert, pic, p, 23.99581946671034)
		})
		test('4', assert => {
			const pic = PICurve.forParametricStartEnd(
				new ProjectedCurveSurface(
					new BezierCurve(
						V(4, 13.8, 15.7),
						V(3, 13.5, 13),
						V(1, 12.8, 6.9),
						V(0, 12.5, 4.200000000000001),
						-0.1,
						1.1,
					),
					V(-1, -0.8, -5.5),
					0,
					1,
					-1,
					0,
				),
				new SemiCylinderSurface(
					new SemiEllipseCurve(
						V(2.5, 13.4, 11.65),
						V(-0.2, -0.16000000000000003, -1.1),
						V(0.2, 0.08000000000000002, 0.68),
						0,
						3.141592653589793,
					),
					V(6.123233995736766e-17, 0.30000000000000004, 2.1000000000000005),
					0,
					3.141592653589793,
					0,
					2,
				),
				V(0.45794084603335994, -0.43492609705736124, 0),
				V(0.5433127907407138, -0.5455452074045032, 0),
				0.05,
				V(-0.0038210833733155064, -0.049853779413943866, 8.673617379884035e-19),
				0,
				11,
			)
			const p = V(2.700000000000252, 13.710000000000084, 13.80000000000074)
			testPointT(assert, pic, p, 6.3461388249955215)
		})
		test('5', assert => {
			const pic = PICurve.forParametricStartEnd(
				new ProjectedCurveSurface(
					new BezierCurve(
						V(0.3333333333333333, -0.1111111111111111, -1.333333333333333),
						V(0.3333333333333333, 0.011618833295731858, -1.333333333333333),
						V(0.2338410555179541, 0.1111111111111111, -1.333333333333333),
						V(0.11111111111111112, 0.1111111111111111, -1.333333333333333),
						0,
						1,
					),
					V(0, 0, -2.222222222222222),
					0,
					1,
					-1,
					0,
				),
				new SemiEllipsoidSurface(
					V3.O,
					V(0.9999999999999999, 0, 0),
					V(0, 0.9999999999999999, 0),
					V(0, 0, -0.9999999999999999),
				),
				V(1.002075784554904, -0.15555163884043577, 0),
				V(-0.047646525196279994, -0.17960385445412552, 2.710505431213761e-20),
				-0.05,
				V(-0.04999136396194676, -0.000929262731508745, 2.710505431213761e-20),
				0,
				21,
			)
			const p = V(0.30360395660901607, -6.873697213268433e-13, -0.9527983194419218)
			testPointT(assert, pic, p, 13.449200584902428)
		})
		test('6', assert => {
			const pic = PICurve.forParametricStartEnd(
				new ProjectedCurveSurface(
					new BezierCurve(
						V(0, 2, 0),
						V(0.3333333333333333, 1, 0),
						V(0.6666666666666666, -1, 0),
						V(1, -2, 0),
						-0.1,
						1.1,
					),
					V(0, 0, -1),
					-0.1,
					1.1,
					-1,
					0,
				),
				new SemiCylinderSurface(
					new SemiEllipseCurve(
						V(0, 0, 0.5),
						V(1.2246467991473533e-17, 0, -0.2),
						V(0, 0.2, 0),
						0,
						3.141592653589793,
					),
					V(1, 0, 6.123233995736766e-17),
					0,
					3.141592653589793,
					0,
					2,
				),
				V(0.5103488594210569, -0.3054969051428064, 0),
				V(0.49441919231892134, -0.6984170426639139, 0),
				0.05,
				V(-0.03401735422540411, 0.03664450315536262, 0),
				0,
				11,
			)
			const p = V(0.49999999999999967, 0, 0.30000000000000004)
			testPointT(assert, pic, p)
		})
	})
	test('isTsWithSurface', assert => {
		const pcs = new ProjectedCurveSurface(BezierCurve.EX2D, V3.Z, undefined, undefined, -2, 2)
			.scale(0.2, 0.2, 1)
			.rotateX(-90 * DEG)
		const ses = SemiEllipsoidSurface.UNIT
		const pic = ses.isCurvesWithSurface(pcs)[0]
		testISTs(assert, pic, new PlaneSurface(P3.XY), 1)
	})
	test('testCurve', assert => {
		const curve = PICurve.forParametricStartEnd(
			new ProjectedCurveSurface(
				new BezierCurve(
					V(0.30000000000000004, -0.1, -1.2),
					V(0.30000000000000004, 0.010456949966158674, -1.2),
					V(0.2104569499661587, 0.1, -1.2),
					V(0.10000000000000002, 0.1, -1.2),
					0,
					1,
				),
				V(0, 0, 2),
				0,
				1,
				0,
				1,
			),
			new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V(0, 0, -1)),
			V(1.0013147208257906, 0.10500326852305793, 0),
			V(-0.04846385366656537, 0.12648109268323876, 0),
			0.05,
			V(-0.049993023321055666, 0.0008352360267516577, 0),
			0,
			13.433240755090269,
		)
		testCurve(assert, curve)
		testCurve(assert, curve.transform(M4.FOO))
	})
	test('testCurve 2', assert => {
		const curve = PICurve.forParametricStartEnd(
			new ProjectedCurveSurface(
				new BezierCurve(
					V(0, 2, 0),
					V(0.3333333333333333, 1, 0),
					V(0.6666666666666666, -1, 0),
					V(1, -2, 0),
					-0.1,
					1.1,
				),
				V(0, 0, -1),
				-0.1,
				1.1,
				-1,
				0,
			),
			new SemiCylinderSurface(
				new SemiEllipseCurve(
					V(0, 0, 0.5),
					V(1.2246467991473533e-17, 0, -0.2),
					V(0, 0.2, 0),
					0,
					3.141592653589793,
				),
				V(1, 0, 6.123233995736766e-17),
				0,
				3.141592653589793,
				0,
				2,
			),
			V(0.5103488594210569, -0.3054969051428064, 0),
			V(0.49441919231892134, -0.6984170426639139, 0),
			0.02,
			V(-0.03401735422540411, 0.03664450315536262, 0),
			2,
			10,
		)
		testCurve(assert, curve, false)
		testCurve(assert, curve.transform(M4.FOO), false)
	})
})
