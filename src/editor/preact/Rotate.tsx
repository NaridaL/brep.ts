///<reference path="../editor.ts"/>
///<reference path="../Sketch.ts"/>
///<reference path="../modes.ts"/>

const {h, render, Component} = preact

const snAndNames = {
	[Extrude.name]: ['EXTR', MODES.EXTRUDE],
	[Rotate.name]: ['ROTA', MODES.ROTATE],
	[Sketch.name]: ['SKTC', MODES.SKETCH],
	[Pattern.name]: ['PTRN', MODES.PATTERN],
	[PlaneDefinition.name]: ['PLNE', MODES.PLANE_DEFINITION],
}

class RotateEditor extends Component<{feature: Rotate, done: any, del: any}, any> {
	render({feature, done, del, notifier}) {
		return h(FeatureEditor, {done, del, notifier, feature},
			<SketchLoopSelect notifier={notifier} feature={feature} propName='_segmentName' label='Loop' />,
			<LineSelect notifier={notifier} feature={feature} propName='axis' label='Axis' />,
			<FeatureAngleEdit notifier={notifier} feature={feature} propName='start' label='Start' />,
			<FeatureAngleEdit notifier={notifier} feature={feature} propName='end' label='End' />,
			<BrepOpSelect notifier={notifier} feature={feature} propName='operation' />)
	}
}

class PatternEditor extends Component<{feature: Pattern, done: any, del: any}, any> {
	render({feature, done, del, notifier}) {
		return h(FeatureEditor, {done, del, notifier, feature},
			<FeaturesSelect notifier={notifier} feature={feature} propName='features' label='Features' />,
			<DirectionSelect notifier={notifier} feature={feature} propName='direction' label='Direction' />,
			<DimensionEdit notifier={notifier} feature={feature} propName='totalLength' label='Total' />,
			<DimensionEdit notifier={notifier} feature={feature} propName='intervalLength' label='Spacing' />,
			<DimensionEdit notifier={notifier} feature={feature} propName='count' label='Count' />)
	}
}

class ExtrudeEditor extends Component<{feature: Extrude, done: any, del: any}, any> {
	render({feature, done, del, notifier}) {
		return h(FeatureEditor, {done, del, notifier, feature},
			<SketchLoopSelect notifier={notifier} feature={feature} propName='_segmentName' label='Loop' />,
			<DimensionEdit notifier={notifier} feature={feature} propName='start' label='Start' />,
			<DimensionEdit notifier={notifier} feature={feature} propName='end' label='End' />,
			<BrepOpSelect notifier={notifier} feature={feature} propName='operation' />)
	}
}

class PlaneEditor extends Component<{feature: Extrude, done: any, del: any}, any> {
	/*
	 <select class="select-select" name="planeType" style="display: inline-block;">
	 <option value="face">Face</option>
	 <option value="points">Through Points</option>
	 <option value="anchorNormal">Anchor &amp; Normal</option>
	 <option value="immediate">Immediate</option>
	 </select>
	 <div><span>Select face</span>
	 <label data-tooltip="Angle blah blub" title="ARFHAKJ"><span>Offset</span>
	 <input class="dimension-input" data-default-value="0" data-feature-property="offset">
	 </label>
	 <label><span>Angle</span>
	 <input class="dimension-input" data-default-value="0" data-feature-property="angle">
	 </label>
	 <label><span>Flip</span>
	 <input class="boolean-input" data-default-value="false" data-feature-property="flipped" type="checkbox">
	 </label>
	 </div>
	 <div>Input
	 <input class="select-text" name="source">
	 </div>
	 </div>
	 */
	render({feature, done, del, notifier}) {
		return h(FeatureEditor, {done, del, notifier, feature},
			<SketchLoopSelect notifier={notifier} feature={feature} propName='_segmentName' />,
			<DimensionEdit notifier={notifier} feature={feature} propName='start' label='Start' />,
			<DimensionEdit notifier={notifier} feature={feature} propName='end' label='End' />,
			<BrepOpSelect notifier={notifier} feature={feature} propName='operation' />)
	}
}

class SketchEditor extends Component<{feature: Extrude, done: any, del: any}, any> {
	render({feature, done, del, notifier}) {
		return h(FeatureEditor, {done, del, notifier, feature},
			<PlaneSelect notifier={notifier} feature={feature} propName='planeRef' />)
	}
}

class FeatureEditor extends Component<any, any> {
	render({done, del, children, notifier, feature, ...props}) {
		const color = '#ff6b6d'
		const style = `display: flex; justify-content: space-between; align-items: center; margin: -2px; background-color: ${color}; padding: 4px;`
		return <div class='editorBox'>

			<div style={style}>
				<strong>{snAndNames[feature.constructor.name][0]}</strong>
				<input type='text' style='width: 70%' value={feature.name} onChange={(e: any) => notifier(feature, 'name', e.target.value)} />
			</div>
			{children}
			<div style='display: flex; justify-content: space-between;'>
				<button style='flex: 1; margin: 2px;' onClick={done}>Done</button>
				<button style='flex: 1; margin: 2px;' onClick={del}>Delete</button>
			</div>
		</div>
	}
}
/*

 div#patternEditor.editorBox
 div(style='display: flex; justify-content: space-between; align-items: center; margin: 2px')
 strong PTRN
 input.string-id-input(type='text' style='width: 70%' data-feature-property='name')
 label
 span.selector.features-select(data-feature-property='features') select features
 div#direction-tooltip(style='display: none;')
 | Select a direction by manually entering a vector or selecting
 ul
 li 1. straight edge: edge direction
 li a curved edge and a point: tangent at point
 li a plane: plane normal
 li 1. non-planar face and 1. point: face normal at point
 div(style='display: flex; justify-content: space-between;' data-tooltipsrc='direction-tooltip' data-tooltipside='left')
 span Direction
 select(data-feature-property='type')
 option(value='auto' selected) Automatic
 option(value='immediate') Immediate
 input.boolean-input(type='checkbox' data-feature-property='directionFlipped')
 div#buh.selector.direction-select(data-feature-property='direction') select things
 div.direction-select-info rewscvhjg
 div.direction-select-things
 div(style='display: flex; justify-content: space-between;align-items: center; margin: 2px;')
 input(type='radio' value='count')
 span(style='flex: 1') count
 input.dimension-input(type='number' data-feature-property='count')
 div(style='display: flex; justify-content: space-between;align-items: center; margin: 2px;')
 input(type='radio' value='intervalLength')
 span(style='flex: 1') spacing
 input.dimension-input(type='number' data-feature-property='intervalLength')
 div(style='display: flex; justify-content: space-between;align-items: center; margin: 2px;')
 input(type='radio' value='totalLength')
 span(style='flex: 1') total length
 input.dimension-input(type='number' data-feature-property='totalLength')
 |
 div(style='display: flex; justify-content: space-between;')
 button(name='done' style='flex: 1; margin: 2px;') Done
 button(name='delete' style='flex: 1; margin: 2px;') Delete
 */

class FeatureAngleEdit extends Component<any, any> {
	static rad2deg(angle: raddd): number {
		const valueInDegrees = angle / DEG
		const rounded = round10(valueInDegrees, -8)
		return (rounded * DEG == angle) ? rounded : valueInDegrees
	}

	onChange = (e: any) => {
		const {feature, notifier, propName} = this.props
		notifier(feature, propName, e.target.value * DEG)
	}

	onWheel = (e: WheelEvent) => {
		const {feature, notifier, propName} = this.props
		const curentValueInDegrees = FeatureAngleEdit.rad2deg(feature[propName])
		const delta = e.shiftKey ? 45 : e.ctrlKey ? 0.1 : 5
		let newValueInDegrees = curentValueInDegrees + -sign(e.deltaY) * delta
		if (round10(curentValueInDegrees, -8) * DEG == feature[propName]) {
			// i.e. current value is rounded, then round this value too
			newValueInDegrees = round10(newValueInDegrees, -8)
		}
		notifier(feature, propName, newValueInDegrees * DEG)
		e.preventDefault()
	}

	render(props) {
		const angle = props.feature[props.propName]
		const valueInDegrees = FeatureAngleEdit.rad2deg(angle)
		return <OneLineThing label={props.label}>
				<input type='text' onChange={this.onChange} value={'' + valueInDegrees} onWheel={this.onWheel}
				       style='width: 50%; text-align: right;' />
			</OneLineThing>
	}
}

class DimensionEdit extends Component<any, any> {
	onChange = (e: any) => {
		const {feature, notifier, propName} = this.props
		notifier(feature, propName, e.target.value)
	}

	onWheel = (e: WheelEvent) => {
		const {feature, notifier, propName} = this.props
		const curentValue = FeatureAngleEdit.rad2deg(feature[propName])
		const delta = e.shiftKey ? 10 : e.ctrlKey ? 0.1 : 1
		let newValue = feature[propName] + -sign(e.deltaY) * delta
		if (round10(curentValue, -8) == feature[propName]) {
			// i.e. current value is rounded, then round this value too
			newValue = round10(newValue, -8)
		}
		notifier(feature, propName, newValue)
		e.preventDefault()
	}

	render(props) {
		const value = props.feature[props.propName]
		//return h('label', {},
		//	props.label,
		//	h('input', {type: 'text', onChange: this.onChange, value: '' + value, onWheel: this.onWheel}))
		return <OneLineThing label={props.label}>
			<input type='text' onChange={this.onChange} value={'' + value} onWheel={this.onWheel}
			style='width: 50%; text-align: right;' />
		</OneLineThing>
	}
}

function OneLineThing(props, state) {
	return <div style='display: flex; justify-content: space-between; align-items: center; margin: 2px'>
		<div style='font-size: small'>{props.label}</div>
		{props.children}
	</div>
}

/*
el.getElements('.face-select, .plane-select, .segment-select, .direction-select, .features-select')
	.removeEvents()
	.removeClass('selecting')
	.addEvent('mouseover', function (e) {
		if (this.linkRef && this.linkRef.get) {
			const target = this.linkRef.get()
			if (target) {
				hoverHighlight = target
			} else {
				missingEls.push(this.linkRef.lastHit)
			}
			paintScreen()
		}
	})
	.addEvent('mouseout', function (e) {
		hoverHighlight = null
		this.linkRef && missingEls.remove(this.linkRef.lastHit)
		paintScreen()
	})

el.getElements('.segment-select')
	.addEvent('click', function (e) {
		const selector = this
		this.addClass('selecting')
		selector.set('text', 'Click on a sketch segment')
		modePush(MODES.SELECT_SEGMENT, function (segmentRef) {
			selector.removeClass('selecting')
			selector.set('text', segmentRef.ref)
			selector.linkRef = segmentRef

			feature[selector.dataset.featureProperty] = segmentRef
			rebuildModel()
		})
	})
*/
class SketchLoopSelect extends Component<any, any> {
	constructor(props) {
		super(props)
		this.state = {selecting: false}
	}
	onClick = e => {
		const {feature, notifier, propName} = this.props
		if (!this.state.selecting) {
			modePush(MODES.SELECT_SEGMENT, segmentRef => {
				notifier(feature, propName, segmentRef)
				modePop()
				this.setState({selecting: false})
			})
		}
		this.setState({selecting: !this.state.selecting})
	}
	render({label, feature, propName}, state) {
		//const selecting = modeStack.some(mode => mode == MODES.SELECT_SEGMENT && )
		const value = feature[propName]
		const classes = classNames('selector', {selecting: state.selecting})
		return <OneLineThing label={label}>
			<div class={classes} onClick={this.onClick} style='display: inline-block; text-align: center; width: 70%'>
				{(state.selecting ? 'Click on a segment' : value ? <RefItem value={value}/> : 'undefined')}
			</div>
		</OneLineThing>
	}
}
class RefItem  extends Component<any, any> {
	onMouseEnter = e => {
		const target = this.props.value.get()
		if (target) {
			hoverHighlight = target
		} else {
			missingEls.push(this.props.value.lastHit)
		}
		paintScreen()
	}
	onMouseLeave = e => {
		hoverHighlight = undefined
		missingEls.remove(this.props.value.lastHit)
		paintScreen()
	}
	render({value, children, ...props}, state) {
		return <div {...props} onMouseLeave={this.onMouseLeave} onMouseEnter={this.onMouseEnter}>{value.ref}</div>
	}
}
class PlaneSelect extends Component<any, any> {
	constructor(props) {
		super(props)
		this.state = {selecting: false}
	}
	onClick = e => {
		const {feature, notifier, propName} = this.props
		if (!this.state.selecting) {
			modePush(MODES.SELECT_PLANE, segmentRef => {
				notifier(feature, propName, segmentRef)
				modePop()
				this.setState({selecting: false})
			})
		}
		this.setState({selecting: !this.state.selecting})
	}
	render(props, state) {
		//const selecting = modeStack.some(mode => mode == MODES.SELECT_SEGMENT && )
		const value = props.feature[props.propName]
		return h('div.selector', {class: classNames({selecting: state.selecting}), onClick: this.onClick},
			'' + (state.selecting ? 'Click on a plane' : value ? value.ref : 'undefined'))
	}
}


abstract class MagicSelect extends Component<any, any> {
	constructor(props) {
		super(props)
		this.state = {selecting: false}
	}

	onClick = e => {
		const {feature, notifier, propName} = this.props
		if (!this.state.selecting) {
			modePush(this.constructor.mode, refs => {
				notifier(feature, propName, refs)
				this.setState({selecting: false})
			})
		}
		this.setState({selecting: !this.state.selecting})
	}
	abstract onMouseOver: (e) => void
	abstract onMouseOut: (e) => void
	render(props, state) {
		//const selecting = modeStack.some(mode => mode == MODES.SELECT_SEGMENT && )
		const refs = props.feature[props.propName]
		const classes = classNames({selecting: state.selecting})
		return <div class={classes} onClick={this.onClick} onMouseOver={this.onMouseOver} onMouseOut={this.onMouseOut}>
			<div style='font-size: small'>{props.label}</div>
			{0 == refs.length
				? 'nothing selected'
				: refs.map(ref => <RefItem value={ref} />)}
		</div>
	}
}
class DirectionSelect extends MagicSelect {
	static mode = MODES.SELECT_DIRECTION
	hit
	onMouseOver = e => {
		const {feature, propName} = this.props
		const whats = feature[propName].map(what => what.get())
		if (whats.every(x => x)) {
			const lineInfo = this.constructor.mode.magic(whats)
			if (lineInfo) {
				const [line, desc] = lineInfo
				missingEls.push(this.hit = line)
			}
		}
	}
	onMouseOut = e => {
		this.hit && missingEls.remove(this.hit)
	}
}
class LineSelect extends MagicSelect {
	static readonly mode = MODES.SELECT_LINE
	hit
	onMouseOver = e => {
		const {feature, propName} = this.props
		const whats = feature[propName].map(what => what.get())
		if (whats.every(x => x)) {
			const lineInfo = this.constructor.mode.magic(whats)
			if (lineInfo) {
				const [line, desc] = lineInfo
				missingEls.push(this.hit = line)
			}
		}
	}
	onMouseOut = e => {
		this.hit && missingEls.remove(this.hit)
	}
}

class FeaturesSelect extends Component<any, any> {
	constructor(props) {
		super(props)
		this.state = {selecting: false}
	}

	onClick = e => {
		const {feature, notifier, propName} = this.props
		if (!this.state.selecting) {
			modePush(MODES.SELECT_FEATURE, refs => {
				notifier(feature, propName, refs)
			})
		}
		this.setState({selecting: !this.state.selecting})
	}

	render({label, ...props}, state) {
		//const selecting = modeStack.some(mode => mode == MODES.SELECT_SEGMENT && )
		const refs = props.feature[props.propName]
		const classes = classNames({selecting: state.selecting})
		return <div class={classes} onClick={this.onClick}
		            data-tooltipside='left'
		            style='display: flex; justify-content: space-between; align-items: top; margin: 2px;'
		            data-tooltip='Select features by clicking on them in the model or the feature stack.'>
			<div>{label}</div>
			<div style='width: 70%; text-align: center'>{0 == refs.length
				? 'nothing selected'
				: refs.map(ref => <RefItem value={ref} />)}
			</div>
		</div>
	}
}

class BrepOpSelect extends Component<any, any> {
	onChange = e => {
		const {feature, notifier, propName} = this.props
		notifier(feature, propName, e.target.value)
	}

	render(props, state) {
		const value = props.feature[props.propName]
		return <select onChange={this.onChange} value={value}>
			<option value='minus'>Difference</option>
			<option value='plus'>Union</option>
			<option value='and'>Intersection</option>
		</select>
	}
}

class FeatureStackDisplay extends Component<any, any> {
	toggleHide(feature) {
		feature.hide = !feature.hide
		paintScreen()
		updateFeatureDisplay()
	}
	featureMouseOver(feature) {
		//const dependencies = featureDependencies(feature)
		//const dependents = featureDependents(feature)
		//div.getChildren().filter((subDiv: any) => dependencies.includes(subDiv.featureLink)).addClass('isDependedOn')
		//div.getChildren().filter((subDiv: any) => dependents.includes(subDiv.featureLink)).addClass('hasDependents')
		hoverHighlight = feature
		updateFeatureDisplay()
		paintScreen()
	}
	onmouseout(e) {
		div.getChildren().removeClass('isDependedOn').removeClass('hasDependents')
	}

	render({features, ...props}, state) {
		const dependencies = hoverHighlight instanceof Feature && featureDependencies(hoverHighlight) || []
		const dependents = hoverHighlight instanceof Feature && featureDependents(hoverHighlight) || []
		const children = []
		for (let featureIndex = 0; featureIndex < features.length; featureIndex++) {
			const feature = features[featureIndex]
			if (rebuildLimit == featureIndex) {
				children.push(<div class='rebuildLimit'>REBUILD LIMIT</div>)
			}
			const classes = classNames({
				'feature': true,
				'isDependentOn': dependencies.includes(feature),
				'hasDependents': dependents.includes(feature)})
			children.push(<div class={classes} onMouseOver={e => this.featureMouseOver(feature)}>
				<span style='width: 42px; display: inline-block;'>{snAndNames[feature.constructor.name][0]}</span>
				<span style='width: 105px; display: inline-block;color: purple;'>{feature.name}</span>

				<a class='icon-link featureDelete' name='delete' onClick={e => featureDelete(feature)}
				   data-tooltip='Delete feature. Dependent features will fail to build and will need to be modified.'>❌</a>

				<a class='icon-link featureEdit' name='edit'
				   onClick={e => modePush(snAndNames[feature.constructor.name][1], feature)}
				   data-tooltip='Modify this feature.'>🖉</a>

				<a class={classNames('icon-link', {'hidden': feature.hide})}
				   onClick={e => this.toggleHide(feature)}
				   data-tooltip='Temporarily prevent this feature from being built. Dependent features will not be built.'>👁</a>

				<a class='icon-link' name='rollBack' onClick={e => featureRollBack(feature, featureIndex)}
				   data-tooltip='Rollback feature.'>R</a>

				{featureError && featureError.feature == feature &&
					<span style='cursor:help; color:red; font-weight: bold;' name='error'
					data-tooltip={featureError.error + featureError.error.stack}>⚠</span>}
			</div>)
		}
		return h('div', {}, children)
	}
}