abstract class Feature {
	abstract apply(publish: (name: string, value: any) => void, color: GL_COLOR, m4?: M4, genFeature?: Feature): void
    name: string
    abstract dependentOnNames(): NameRef[]
}