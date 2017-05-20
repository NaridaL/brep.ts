class PatternEditor extends React.Component<Pattern & {magic: any}, any> {
	constructor() {
		super()
	}

	render() {
		return (<div id="patternEditor" className='editorBox'>
			<span className='header'>PTRN</span>
			<input type='text' data-feature-property='name' className='string-id-input' style={{display: 'inline-block', width: '100px'}} />
			<label>
				<span>Feature</span>
				<span className='feature-select selector' data-feature-property='feature'>sel feat</span></label>

			{this.props.dimensions.map(dim => React.createElement(PatternDimension, {...dim, magic: this.props.magic}))}

			<button onClick={e => (this.props.magic(x => x.addDimension()), this.forceUpdate())}>+ Dimension</button>

			<button name='done'>Done</button>
			<button name='delete'>Delete</button>
		</div>)
	}
}

function PatternDimension(props: PatternDimensionType) {
	return <div>
		<label>
			<span>Direction</span>
			<span className='direction-select selector' data-feature-property='direction'>sel dir</span></label>
		<label>
			<span>Flip</span>
			<input className='boolean-input' data-default-value='false' data-feature-property='directionFlipped'
			       type='checkbox' /></label>

		{[['count', 'count', 2], ['totalLength', 'total length', 20], ['intervalLength', 'spacing', 10]]
			.map(([featProp, name]) =>
				<div>
					<input type='radio' name='test_0' style={{display: 'inline-block'}}
					       onClick={e => (this.props.magic(x => x.dimensions[0].passive = featProp), this.forceUpdate())}/>
					<label style={{display: 'inline-block'}}>
						<span>{name}</span>
						<input type='text' className='dimension-input' disabled={props.passive == featProp}
						       data-feature-property={featProp} data-default-value='0'
						       style={{display: 'inline-block', width: '40px'}} /></label></div>)}</div>
}

