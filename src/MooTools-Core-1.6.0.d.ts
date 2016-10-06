interface ElementConstructor {
	prototype: Element
	new(): Element
	new (tagNameOrCSSSelector: string, properties: {[prop: string]: any})
}

/////// Element
interface Element {
	new (tagNameOrCSSSelector: string, properties: {[prop: string]: any})
	new (el: HTMLElement, properties: {})

	/**
	 * Gets the first descendant element whose tag name matches the tag provided. CSS selectors may also be passed.
	 * @param tag Tag name of the element to find or a CSS Selector.
	 */
	getElement<T extends HTMLElement>(tag: string): T

	/**
	 * Collects all descendant elements whose tag name matches the tag provided. CSS selectors may also be passed.
	 * @param tag Tag name of the element to find or a CSS Selector.
	 */
	getElements(tag: string): Elements

	/**
	 * Gets the element with the specified id found inside the current Element.
	 *
	 * Notes: This method is not provided for Document instances as document.getElementById is provided natively.
	 *
	 * @param id The ID of the Element to find.
	 * @returns If a match is found, returns that Element. Otherwise, returns null.
	 */
	getElementById(id: string)

	/**
	 * This is a "dynamic arguments" method. Properties passed in can be any of the 'set' properties in the
	 * Element.Properties Object.
	 *
	 *
	 Two Arguments (property, value)
	 property:string The string key from the Element.Properties Object representing the property to set.
	 value:mixed The value to set for the specified property.
	 One Argument (properties)
	 properties:object Object with its keys/value pairs representing the properties and values to set for the Element (as described below).
	 Returns:
	 (element) This Element.
	 Examples:
	 With Property and Value:
	 $('myElement').set('text', 'text goes here');
	 $('myElement').set('class', 'active');
	 // the 'styles' property passes the object to Element:setStyles.
	 var body = $(document.body).set('styles', {
	font: '12px Arial',
	color: 'blue'
});
	 With an Object:
	 var myElement = $('myElement').set({
		// the 'styles' property passes the object to Element:setStyles.
		styles: {
			font: '12px Arial',
			color: 'blue',
			border: '1px solid #f00'
		},
		// the 'events' property passes the object to Element:addEvents.
		events: {
			click: function(){ alert('click'); },
			mouseover: function(){ this.addClass('over'); }
		},
		//Any other property uses Element:setProperty.
		id: 'documentBody'
	});
	 Notes:
	 All the property arguments are passed to the corresponding method of the Element.Properties Object.
	 If no matching property is found in Element.Properties, it falls back to Element:setProperty.
	 Whenever using Element:setProperty to set an attribute, pass in the lowercase, simplified form of the property. For example:
	 use 'for', not 'htmlFor',
	 use 'class', not 'className'
	 use 'frameborder', not 'frameBorder'
	 etc.
	 In IE8 or lower, it is not possible to set type multiple times. It will throw an error.
	 See Also:
	 Element, Element.Properties, Element:setProperty, Element:addEvents, Element:setStyles
	 */
	set(property: string, value: any)
	set(properties: { [property: string] : any })

	/**
	 *
	 Element Method: get
	 Back to Top
	 This is a "dynamic arguments" method. Properties passed in can be any of the 'get' properties in the Element.Properties Object.
	 Syntax:
	 myElement.get(property);
	 Arguments:
	 property:string The string key from the Element.Properties Object representing the property to get.
	 Returns:
	 (mixed) The result of calling the corresponding 'get' function in the Element.Properties Object.
	 Examples:
	 Using Custom Getters:
	 var tag = $('myDiv').get('tag'); // returns "div".
	 Fallback to Element Attributes:
	 var id = $('myDiv').get('id'); // returns "myDiv".
	 var value = $('myInput').get('value'); // returns the myInput element's value.
	 Notes:
	 If the corresponding accessor doesn't exist in the Element.Properties Object, the result of Element:getProperty on the property passed in is returned.
	 See Also:
	 Element, Element.Properties, Element:getProperty
	 */

	/**
	 This is a "dynamic arguments" method. Properties passed in can be any of the 'erase' properties in the Element.Properties Object.
	 Syntax:
	 myElement.erase(property);
	 Arguments:
	 property:string The string key from the Element.Properties Object representing the property to erase.
	 Returns:
	 (mixed) The result of calling the corresponding 'erase' function in the Element.Properties Object.
	 Examples:
	 $('myDiv').erase('id'); //Removes the id from myDiv.
	 $('myDiv').erase('class'); //myDiv element no longer has any class names set.
	 Note:
	 If the corresponding eraser doesn't exist in the Element.Properties Object, Element:removeProperty is called with the property passed in.
	 See Also:
	 Element, Element.Properties, Element:removeProperty
	 */
	erase(property: string)

	/**
	 *
	 * @returns If the element matched, returns true. Otherwise, returns false.
	 */
	match(match: string): boolean
	match(match: Element): boolean

	contains(el: Element): boolean

	inject(elOrId: Element | string, where?: 'top' | 'bottom' | 'after' | 'before'): this
	grab(elOrId: Element | string, where?: 'top' | 'bottom' | 'after' | 'before'): this

	adopt(elsOrId: Element | Element[] | string, ...others: (Element | Element[])[]): this

	wraps(elOrId: Element | string, where: 'top'|'bottom'): this

	addEvent(type: string, fn: (e: DOMEvent) => any):this
	addEvents(events:{[type: string] :  (e: DOMEvent) => any}):this

	removeEvents(type:string):this
	removeEvents(events:{[type: string] :  (e: DOMEvent) => any}):this

	/**
	 Works like [Element:grab](#Element:grab), but instead of accepting an id or an element, it only accepts an HTML string.
	 The HTML string will be parsed to create new DOM elements, and then injected relative to the element from where the method
	 was called.
	 ### Examples:
	 ##### HTML
	 <div id="myElement">Hey.</div>
	 ##### JavaScript
	 $('myElement').appendHTML(' <strong>Howdy.</strong>');
	 ##### Resulting HTML
	 <div id="myElement">Hey. <strong>Howdy.</strong></div>
	 ### Notes:
	 - This method does *not* use the `innerHTML` property of an element but instead creates elements
	 directly before injecting them. Thus, it is safe to use in cases where you don't want to destroy
	 any descendant elements already present in the parent.
	 - This method uses `insertAdjacentHTML` when available.
	 ### See Also:
	 - [MDN Element:insertAdjacentHTML][].
	 */
	/**
	 *
	 * @param html The HTML string to append.
	 * @param where The position to inject the text to.  Defaults to 'bottom'.
	 */
	appendHTML(html: string, where?: 'top' | 'bottom' | 'after' | 'before')
	/**

	 Works like [Element:grab](#Element:grab), but instead of accepting an id or an element, it only accepts text.
	 A text node will be created inside this Element, in either the top or bottom position.

	 ### Returns:

	 * (*element*) The current Element instance.

	 ### Examples:

	 ##### HTML

	 <div id="myElement">Hey.</div>

	 ##### JavaScript

	 $('myElement').appendText(' Howdy.');

	 ##### Resulting HTML

	 <div id="myElement">Hey. Howdy.</div>

	 */
	/**
	 *
	 * @param html The text to append.
	 * @param where The position to inject the text to. Defaults to 'bottom'.
	 */
	appendText(text: string, where?: 'top' | 'bottom' | 'after' | 'before')
	/**
	 Removes the Element from the DOM.
	 ### Examples:
	 ##### HTML
	 <div id="myElement"></div>
	 <div id="mySecondElement"></div>
	 ##### JavaScript
	 $('myElement').dispose();
	 ##### Resulting HTML
	 <div id="mySecondElement"></div>
	 ### See Also:
	 - [MDN Element:removeChild][]
	 @returns This Element. Useful to always grab the return from this function, as the element could be
	      [injected](#Element:inject) back.
	 */
	dispose(): this

	/**
	 Clones the Element and returns the cloned one.
	 ### Arguments:
	 @param contents When set to false the Element's contents are not cloned. default=true
	 @param keepid When true the cloned Element keeps the id attribute, if present. Same goes for any of the cloned
	      childNodes. default=false
	 ### Examples:
	 ##### HTML
	 <div id="myElement">ciao</div>
	 ##### JavaScript
	 // clones the Element and appends the clone after the Element.
	 var clone = $('myElement').clone().inject('myElement','after');
	 ##### Resulting HTML
	 <div id="myElement">ciao</div>
	 <div>ciao</div>
	 ### Note:
	 - The returned Element does not have attached events. To clone the events use [Element:cloneEvents](/core/Element/Element.Event#Element:cloneEvents).
	 - Values stored in Element.Storage are not cloned.
	 - The clone element and its children are stripped of ids, unless otherwise specified by the keepid parameter.
	 ### See Also:
	 - [Element:cloneEvents](/core/Element/Element.Event#Element:cloneEvents).
	 */
	clone(contents: boolean, keepid: boolean): Element

	/**
	 Replaces the passed Element with Element.

	 ### Arguments:

	 @param el A string id representing the Element to be replaced, or an Element reference.

	 ### Examples:

	 $('myNewElement').replaces($('myOldElement'));
	 //$('myOldElement') is gone, and $('myNewElement') is in its place.

	 ### See Also:

	 - [MDN Element:replaceChild][]

	 */
	replace(elOrId: string | Element): this

	/**
	 Tests the Element to see if it has the passed in className.

	 ### Arguments:

	 @param className The class name to test.

	 ### Returns:

	 * (*boolean*) Returns true if the Element has the class, otherwise false.

	 ### Examples:

	 ##### HTML

	 <div id="myElement" class="testClass"></div>

	 ##### JavaScript

	 $('myElement').hasClass('testClass'); // returns true

	 */
	/**
	 *
	 * @param className The class name to test.
	 */
	hasClass(className: string): boolean

	/**
	 Adds the passed in class to the Element, if the Element doesn't already have it.

	 ### Arguments:

	 @param className The class name to add.

	 ### Examples:

	 ##### HTML

	 <div id="myElement" class="testClass"></div>

	 ##### JavaScript

	 $('myElement').addClass('newClass');

	 ##### Resulting HTML

	 <div id="myElement" class="testClass newClass"></div>

	 */
	addClass(className: string): this

	/**
	 Works like [Element:addClass](#Element:addClass), but removes the class from the Element.

	 ### Arguments:

	 @param className The class name to remove.

	 ### Examples:

	 ##### HTML

	 <div id="myElement" class="testClass newClass"></div>

	 ##### JavaScript

	 $('myElement').removeClass('newClass');

	 ##### Resulting HTML

	 <div id="myElement" class="testClass"></div>

	 */
	removeClass(className: string): this

	/**
	 Adds or removes the passed in class name to the Element, depending on whether or not it's already present.

	 ### Arguments:

	 @param className The class to add or remove.
	 @param force Force the class to be either added or removed

	 ### Examples:

	 ##### HTML

	 <div id="myElement" class="myClass"></div>

	 ##### JavaScript

	 $('myElement').toggleClass('myClass');

	 ##### Resulting HTML

	 <div id="myElement" class=""></div>

	 ##### JavaScript

	 $('myElement').toggleClass('myClass');

	 ##### Resulting HTML

	 <div id="myElement" class="myClass"></div>

	 */
	toggleClass(className: string, force?: boolean): this

	/**
	 * Returns the previousSibling of the Element (excluding text nodes).
	 *
	 * @param match A comma separated list of tag names to match the found element(s) with. A full CSS selector can be
	 *     passed.
	 * @returns The previous sibling Element or null if none found.
	 */
	getPrevious(match?: string): Element | null

	/**
	 * Like [Element:getPrevious][], but returns a collection of all the matched previousSiblings.
	 *
	 * @param match
	 */
	getAllPrevious(match?: string): Elements

	getNext(match?: string): Element | null

	/**
	 * Like [Element:getPrevious][], but returns a collection of all the matched nextSiblings.
	 * @param match
	 */
	getAllNext(match?: string): Elements

	/**
	 Gets the first element that matches the passed in expression.

	 ### Arguments:

	 1. match - (*string*, optional): A full CSS selector to match the found element(s) with.

	 ### Returns:

	 * (*mixed*) The first found element or null if none found.

	 */
	getFirst(match?: string): Element | null

	/**
	 Gets the last element that matches the passed in expression.

	 ### Arguments:

	 1. match - (*string*, optional): A full CSS selector to match the found element(s) with.

	 ### Returns:

	 * (*mixed*) The last found element, or returns null if none found.

	 */
	getLast(match?: string): Element | null

	/**
	 Works as [Element:getPrevious][], but tries to find the parentNode.

	 ### Arguments:

	 1. match - (*string*, optional): A tag name to match the found element(s) with. A full CSS selector can be passed.

	 ### Returns:

	 * (*mixed*) The target Element's parent or null if no matching parent is found.

	 */
	getParent(match?: string): Element | null

	/**
	 Like [Element:getParent](#Element:getParent), but returns a collection of all the matched parentNodes up the tree.

	 ### Returns:

	 * (*array*) If no matching parents are found, an empty array is returned.

	 */
	getParents(match?: string): Element[]

	/**
	 Like [Element:getAllPrevious][] but returns all Element's previous and next siblings (excluding text nodes). Returns as [Elements][].
	 ### Arguments:
	 1. match - (*string*, optional): A tag name to match the found element(s) with. A full CSS selector can be passed.
	 ### Returns:
	 * (*array*) A [Elements](#Elements) array with all of the Element's siblings, except the text nodes.
	 */
	getSiblings(match?: string): Elements

	/**
	 Returns all the Element's children (excluding text nodes). Returns as [Elements][].
	 ### Arguments:
	 1. match - (*string*, optional): A tag name to match the found element(s) with. A full CSS selector can be passed.
	 ### Returns:
	 * (*array*) A [Elements](#Elements) array with all of the Element's children, except the text nodes.
	 ### Note:
	 The difference between the methods *getChildren* and *getElements* is that getChildren will only return its direct children while getElements searches for all the Elements in any depth.
	 */
	getChildren(match?: string): Elements

	/**
	 Empties an Element of all its children.

	 ### Examples:

	 ##### HTML

	 <div id="myElement">
	 <p></p>
	 <span></span>
	 </div>

	 ##### JavaScript

	 $('myElement').empty();

	 ##### Resulting HTML

	 <div id="myElement"></div>

	 ### Note:

	 This method does not garbage collect the children. Use [Element:destroy][] instead.

	 */
	empty(): this

	/**
	 * Removes the Element and its children from the DOM and prepares them for garbage collection.
	 */
	destroy(): null

	/**

	 Reads the child inputs of the Element and generates a query string based on their values.

	 ### Returns:

	 * (*string*)

	 ### Examples:

	 ##### HTML

	 <form id="myForm" action="submit.php">
	 <input name="email" value="bob@bob.com" />
	 <input name="zipCode" value="90210" />
	 </form>

	 ##### JavaScript

	 $('myForm').toQueryString(); // returns "email=bob@bob.com&zipCode=90210".
	 @returns A string representation of a all the input Elements' names and values.
	 */
	toQueryString(): string

	/**
	 Returns the selected options of a select element.

	 ### Returns:

	 * (*array*) An array of the selected elements.

	 ### Examples:

	 ##### HTML

	 <select id="country-select" name="country">
	 <option value="US">United States</option
	 <option value ="IT">Italy</option>
	 </select>

	 ##### JavaScript

	 $('country-select').getSelected(); // returns whatever the user selected.

	 ### Note:

	 This method returns an array, regardless of the multiple attribute of the select element.
	 If the select is single, it will return an array with only one item.

	 */
	getSelected(): Elements

	/**
	 Returns a single element attribute.

	 ### Arguments:

	 * property - (*string*) The property to be retrieved.

	 ### Returns:

	 * (*string*) A string containing the Element's requested property.

	 ### Examples:

	 ##### HTML

	 <img id="myImage" src="mootools.png" title="MooTools, the compact JavaScript framework" alt="" />

	 ##### JavaScript

	 var imgProps = $('myImage').getProperty('src'); // returns: 'mootools.png'.

	 */
	getProperty(property: string): string

	/**
	 Gets multiple element attributes.

	 ### Arguments:

	 * properties - (*strings*) Any number of properties to be retrieved.

	 ### Returns:

	 * (*object*) An object containing all of the Element's requested properties.

	 ### Examples:

	 ##### HTML

	 <img id="myImage" src="mootools.png" title="MooTools, the compact JavaScript framework" alt="" />

	 ##### JavaScript

	 var imgProps = $('myImage').getProperties('id', 'src', 'title', 'alt');
	 // returns: { id: 'myImage', src: 'mootools.png', title: 'MooTools, the compact JavaScript framework', alt: '' }

	 */
	getProperties(...properties: string[]): { [property: string]: string }

	/**
	 Sets an attribute or special property for this Element.
	 ### Arguments:
	 @param property The property to assign the value passed in.
	 @param value The value to assign to the property passed in.
	 ### Examples:
	 ##### HTML
	 <img id="myImage" />
	 ##### JavaScript
	 $('myImage').setProperty('src', 'mootools.png');
	 ##### Resulting HTML
	 <img id="myImage" src="mootools.png" />
	 ### Note
	 - Whenever using [Element:setProperty][] to set an attribute, pass in the lowercase, simplified form of the property. For example:
	 - use 'for', not 'htmlFor',
	 - use 'class', not 'className'
	 - use 'frameborder', not 'frameBorder'
	 - etc.
	 - When setting the `src` property for an image file, be sure to remove the `width` and `height` attribute (use `Element.removeAttribute`). IE7, and less, set and freeze the `width` and `height` of an image if previously specified.
	 */
	setProperty(property: string, value: any): this

	/**
	 Sets numerous attributes for the Element.

	 ### Arguments:

	 @param properties An object with key/value pairs.

	 ### Examples:

	 ##### HTML

	 <img id="myImage" />

	 ##### JavaScript

	 $('myImage').setProperties({
	 src: 'whatever.gif',
	 alt: 'whatever dude'
	 });

	 ##### Resulting HTML

	 <img id="myImage" src="whatever.gif" alt="whatever dude" />

	 */
	setProperties(properties: { [property: string]: string }): this

	/**
	 Removes an attribute from the Element.

	 ### Arguments:

	 @param property The attribute to remove.

	 ### Examples:

	 ##### HTML

	 <a id="myAnchor" href="#" onmousedown="alert('click');"></a>

	 ##### JavaScript

	 //Eww... inline JavaScript is bad! Let's get rid of it.
	 $('myAnchor').removeProperty('onmousedown');

	 ##### Resulting HTML

	 <a id="myAnchor" href="#"></a>

	 */
	removeProperty(property: string): this

	/**
	 Removes numerous attributes from the Element.

	 ### Arguments:

	 @param properties The attributes to remove, separated by comma.

	 ### Examples:

	 ##### HTML

	 <a id="myAnchor" href="#" title="hello world"></a>

	 ##### JavaScript

	 $('myAnchor').removeProperties('id', 'href', 'title');

	 ##### Resulting HTML

	 <a></a>
	 */
	removeProperties(...properties: string[]): this
	/**
	 Stores an item in the Elements Storage, linked to this Element.

	 ### Arguments:

	 @param key The key you want to assign to the stored value.
	 @param value Any value you want to store.

	 ### Example:

	 $('element').store('someProperty', someValue);

	 */
	store(key: string, value: any): this
	/**
	 Retrieves a value from the Elements storage.

	 ### Arguments:

	 @param key The key you want to retrieve from the storage.
	 @param default Default value to store and return if no value is stored.

	 ### Returns:

	 * (*mixed*) The value linked to the key.

	 ### Example:

	 $('element').retrieve('someProperty'); // returns someValue (see example above)

	 */
	retrieve(key: string): any
	/**
	 Eliminates a key from the Elements storage.

	 ### Arguments:

	 @param key The key you want to eliminate from the storage.

	 ### Example:

	 $('element').eliminate('someProperty');

	 */
	eliminate(key: string): this
}

/////// Element.Styles
interface Element {
	/**
	 *
	 Sets a CSS property to the Element.
	 ### Example:
	 //Both lines have the same effect.
	 $('myElement').setStyle('width', '300px'); // the width is now 300px.
	 $('myElement').setStyle('width', 300); // the width is now 300px.
	 ### Notes:
	 - All number values will automatically be rounded to the nearest whole number.
	 * @param property The property to set.
	 * @param value The value to which to set it. Numeric values of properties requiring a unit will automatically be
	 *     appended with 'px'.
	 */
	setStyle(property:CSSProperty, value:string | number): this
	/**
	 * Returns the style of the Element given the property passed in.
	 *
	 * @example $('myElement').getStyle('width'); // returns "300px".
	 * @example $('myElement').getStyle('width').toInt(); // returns 300.
	 * @param property The css style property you want to retrieve.
	 * @returns
	 * The style value.
	 */
	getStyle(property:CSSProperty): string
	/*
	 Applies a collection of styles to the Element.

	 ### Arguments:

	 1. styles - (*object*) An object of property/value pairs for all the styles to apply.

	 ### Example:

	 $('myElement').setStyles({
	 border: '1px solid #000',
	 width: 300,
	 height: 400
	 });

	 ### See Also:

	 - [Element:getStyle][]


	 */
	setStyles(styles: {[property: 'A' | 'b']: string|number}): this
	/**
	 *
	 ### Examples:

	 $('myElement').getStyles('width', 'height', 'padding');
	 // returns {width: '10px', height: '10px', padding: '10px 0px 10px 0px'}

	 ### See Also:

	 - [Element:getStyle][]
	 * Returns an object of styles of the Element for each argument passed in.
	 * @param properties Any number of style properties.
	 * @returns
	 * An key/value object with the CSS styles as computed by the browser.
	 */
	setStyles(...properties: CSSProperty[]): {[property: string]: string|number}
}

/////// Element.Position
interface Element {
	/**
	 * Type: Element {#Element}
	 =========================
	 Custom Type to allow all of its methods to be used with any DOM element via the dollar function [$][].
	 ### Notes:
	 * These methods don't take into consideration the body element margins and borders. If you need margin/borders on
	 * the body, consider adding a wrapper div, but always reset the margin and borders of body to 0. If you need to
	 * measure the properties of elements that are not displayed (either their display style is none or one of their
	 * parents display style is none), you will need to use [Element.measure][] to expose it.
	 ### Credits:
	 - Element positioning based on the [qooxdoo](http://qooxdoo.org/) code and smart browser fixes, [LGPL License](http://www.gnu.org/licenses/lgpl.html).
	 - Viewport dimensions based on [YUI](http://developer.yahoo.com/yui/) code, [BSD License](http://developer.yahoo.com/yui/license.html).
	 Element Method: scrollTo {#Element:scrollTo}
	 --------------------------------------------
	 Scrolls the element to the specified coordinated (if the element has an overflow).
	 The following method is also available on the Window object.
	 ### Syntax:
	 myElement.scrollTo(x, y);
	 ### Arguments:
	 1. x - (*number*) The x coordinate.
	 2. y - (*number*) The y coordinate.
	 ### Example:
	 $('myElement').scrollTo(0, 100);
	 ### See Also:
	 - [MDN Element:scrollLeft][], [MDN Element:scrollTop][]
	 */
	scrollTo(x: number, y: number): this
	/*


	 Element Method: getSize {#Element:getSize}
	 ------------------------------------------

	 Returns the height and width of the Element, taking into account borders and padding.
	 The following method is also available on the Window object.

	 ### Syntax:

	 myElement.getSize();

	 ### Returns:

	 * (*object*) An object containing the width (as x) and the height (as y) of the target Element.

	 ### Example:

	 var size = myElement.getSize();
	 alert('The element is ' + size.x + ' pixels wide and ' + size.y + 'pixels high.');

	 ### Note:

	 If you need to measure the properties of elements that are not displayed (either their display style is none or one of their parents display style is none), you will need to use [Element.measure][] to expose it.

	 */
	getSize(): {x: number, y: number}
	/*

	 Element Method: getScrollSize {#Element:getScrollSize}
	 ------------------------------------------------------

	 Returns an Object representing the size of the target Element, including scrollable area.
	 The following method is also available on the Window object.

	 ### Syntax:

	 myElement.getScrollSize();

	 ### Returns:

	 * (*object*) An object containing the x and y dimensions of the target Element.

	 ### Example:

	 var scroll = $('myElement').getScrollSize();
	 alert('My element can scroll to ' + scroll.y + 'px'); // alerts 'My element can scroll down to 820px'

	 ### See Also:

	 - [MDN Element:scrollLeft][], [MDN Element:scrollTop][], [MDN Element:offsetWidth][], [MDN Element:offsetHeight][], [MDN Element:scrollWidth][], [MDN Element:scrollHeight][]

	 ### Note:

	 If you need to measure the properties of elements that are not displayed (either their display style is none or one of their parents display style is none), you will need to use [Element.measure][] to expose it.

	 */
	getScrollSize(): {x: number, y: number}
	/*

	 Element Method: getScroll {#Element:getScroll}
	 ----------------------------------------------

	 Returns an Object representing how far the target Element is scrolled in either direction.
	 The following method is also available on the Window object.

	 ### Syntax:

	 myElement.getScroll();

	 ### Returns:

	 * (*object*) An object containing the x and y dimensions of the target Element's scroll.

	 ### Example:

	 var scroll = $('myElement').getScroll();
	 alert('My element is scrolled down ' + scroll.y + 'px'); // alerts 'My element is scrolled down to 620px'

	 ### Note:

	 If you need to measure the properties of elements that are not displayed (either their display style is none or one of their parents display style is none), you will need to use [Element.measure][] to expose it.

	 */
	getScroll(): {x: number, y: number}
	/*


	 Element Method: getPosition {#Element:getPosition}
	 --------------------------------------------------

	 Returns the real offsets of the element.

	 ### Syntax:

	 myElement.getPosition(relative);

	 ### Arguments:

	 relative - (Element, defaults to the document) If set, the position will be relative to this Element.

	 ### Returns:

	 * (*object*) An object with the x and y coordinates of the Element's position.

	 ### Example:

	 $('element').getPosition(); // returns {x: 100, y: 500};

	 ### See Also:

	 - [QuirksMode: Find position](http://www.quirksmode.org/js/findpos.html)

	 ### Note:

	 If you need to measure the properties of elements that are not displayed (either their display style is none or one of their parents display style is none), you will need to use [Element.measure][] to expose it.
	 */
	getPosition(relative?: Element): {x: number, y: number}
	/*

	 Element Method: setPosition {#Element:setPosition}
	 --------------------------------------------------

	 Sets the position of the element's *left* and *top* values to the x/y positions you specify.

	 ### Syntax

	 myElement.setPosition(positions);

	 ### Arguments

	 1. positions - (*object*) an object with x/y values (integers or strings, i.e. 10 or "10px")

	 ### Returns

	 * (*element*) the element that is positioned.

	 ### Example

	 myElement.setPosition({x: 10, y: 100});

	 */
	setPosition(positions: {x: number | string, y: number | string}): this
	/*


	 Element Method: getCoordinates {#Element:getCoordinates}
	 --------------------------------------------------------

	 Returns an object with width, height, left, right, top, and bottom coordinate values of the Element.

	 ### Syntax:

	 myElement.getCoordinates(relative);

	 ### Arguments:

	 relative - (*element*, optional) if set, the position will be relative to this element, otherwise relative to the document.

	 ### Returns:

	 * (*object*) An object containing the Element's current: top, left, width, height, right, and bottom.

	 ### Example:

	 var myValues = $('myElement').getCoordinates();

	 #### Returns:

	 {
		 top: 50,
		 left: 100,
		 width: 200,
		 height: 300,
		 right: 300,
		 bottom: 350
	 }

	 ### See Also:

	 [Element:getPosition](#Element:getPosition)

	 ### Note:

	 If you need to measure the properties of elements that are not displayed (either their display style is none or one of their parents display style is none), you will need to use [Element.measure][] to expose it.
	 */
	getCoordinates(relative?: Element): {
		top: number, left: number, width: number, height: number, right: number, bottom: number}
	/*

	 Element Method: getOffsetParent {#Element:getOffsetParent}
	 ----------------------------------------------------------

	 Returns the parent of the element that is positioned, if there is one.

	 ### Syntax

	 myElement.getOffsetParent();

	 ### Returns

	 * (*mixed*) If the element has a parent that is positioned, it returns that element, otherwise it returns *null*.



	 [$]: /core/Element/Element#Window:dollar
	 [MDN Element:scrollLeft]: https://developer.mozilla.org/en/DOM/element.scrollLeft
	 [MDN Element:scrollTop]: https://developer.mozilla.org/en/DOM/element.scrollTop
	 [MDN Element:offsetWidth]: https://developer.mozilla.org/en/DOM/element.offsetWidth
	 [MDN Element:offsetHeight]: https://developer.mozilla.org/en/DOM/element.offsetHeight
	 [MDN Element:scrollWidth]: https://developer.mozilla.org/en/DOM/element.scrollWidth
	 [MDN Element:scrollHeight]: https://developer.mozilla.org/en/DOM/element.scrollHeight
	 [Element.measure]: /more/Element/Element.Measure

	 */
	getOffsetParent(): Element | null
}

interface Elements extends Element, Array<Element> {
	each(f: (el:Element)=>void):this
    filter(f: (el:Element)=>boolean):Elements
}

interface DOMEvent {

	/**
	 * Original event.
	 */
	event: Event

	/**
	 * The position of the mouse, relative to the full window.
	 */
	page: {x: number, y: number}

	/**
	 * The position of the mouse, relative to the viewport.
	 */
	client: {x: number, y: number}

	/**
	 * True if the user clicked the right mousebutton
	 */
	rightClick: boolean

	/**
	 * The amount of third button scrolling.
	 */
	wheel: number

	/**
	 * The event related target.
	 */
	relatedTarget: Element

	/**
	 * The event target.
	 */
	target: Element

	/**
	 * The keycode of the key pressed.
	 */
	code: number

	/**
	 * The key pressed as a lowercase string. key can be 'enter', 'up', 'down', 'left', 'right', 'space', 'backspace',
	 * 'tab', 'delete', and 'esc'.
	 */
	key: string

	/**
	 * True if the user pressed the shift key.
	 */
	shift: boolean

	/**
	 * True if the user pressed the control key.
	 */
	control: boolean

	/**
	 * True if the user pressed the alt key.
	 */
	alt: boolean

	/**
	 * True if the user pressed the meta key.
	 */
	meta: boolean

	/**
	 * Stop an event from propagating and also executes preventDefault.
	 */
	stop(): this

	/**
	 * Cross browser method to stop the propagation of an event (this stops the event from bubbling up through the DOM).
	 */
	stopPropagation(): this

	/**
	 * Cross browser method to prevent the default action of the event.
	 */
	preventDefault(): this

	/**
	 * This function allows to add an additional event key code.
	 * @param keyCode
	 * @param keyName
	 */
	defineKey(keyCode: number, keyName: string)

}
declare function $$(selector: string): Elements

/**
 *
 Function: document.id {#Window:document-id}
 -------------------------------------------
 The document.id function has a dual purpose: Getting the element by its id, and making an element in Internet Explorer "grab" all the [Element][] methods.
 ### Syntax:
 var myElement = document.id(el);
 ### Arguments:
 1. el - The Element to be extended. Can be one of the following types:
 * (*element*) The element will be extended if it is not already.
 * (*string*) A string containing the id of the DOM element desired.
 * (*object*) If the object has a toElement method, toElement will be called to get the Element.
 ### Returns:
 * (*element*) A DOM element.
 * (*null*) Null if no matching id was found or if toElement did not return an element.
 ### Examples:
 #### Get a DOM Element by ID:
 var myElement = document.id('myElement');
 #### Get a DOM Element by reference:
 var div = document.getElementById('myElement');
 div = document.id(div); // the element with all the Element methods applied.
 ### Notes:
 - This method is useful when it's unclear if working with an actual element or an id.  It also serves as a shorthand for document.getElementById().
 - In Internet Explorer, the [Element][] is extended the first time document.id is called on it, and all the [Element][] Methods become available.
 - Browsers with native HTMLElement support, such as Safari, Firefox, and Opera, apply all the [Element][] Methods to every DOM element automatically.
 - Because MooTools detects if an element needs to be extended or not, this function may be called on the same Element many times with no ill effects.
 */
declare function $<T extends Element>(selector: string): T | null
declare function $(elementable: { toElement: () => HTMLElement }): Element | null
declare function $(anything:any|null): null
declare function $<T extends Element>(element: T): T



type CSSProperty = 'azimuth'
	| 'backgroundAttachment'
	| 'backgroundColor'
	| 'backgroundImage'
	| 'backgroundPosition'
	| 'backgroundRepeat'
	| 'background'
	| 'borderCollapse'
	| 'borderColor'
	| 'borderSpacing'
	| 'borderStyle'
	| 'borderTop'
	| 'borderRight'
	| 'borderBottom'
	| 'borderLeft'
	| 'borderTopColor'
	| 'borderRightColor'
	| 'borderBottomColor'
	| 'borderLeftColor'
	| 'borderTopStyle'
	| 'borderRightStyle'
	| 'borderBottomStyle'
	| 'borderLeftStyle'
	| 'borderTopWidth'
	| 'borderRightWidth'
	| 'borderBottomWidth'
	| 'borderLeftWidth'
	| 'borderWidth'
	| 'border'
	| 'bottom'
	| 'captionSide'
	| 'clear'
	| 'clip'
	| 'color'
	| 'content'
	| 'counterIncrement'
	| 'counterReset'
	| 'cueAfter'
	| 'cueBefore'
	| 'cue'
	| 'cursor'
	| 'direction'
	| 'display'
	| 'elevation'
	| 'emptyCells'
	| 'float'
	| 'fontFamily'
	| 'fontSize'
	| 'fontStyle'
	| 'fontVariant'
	| 'fontWeight'
	| 'font'
	| 'height'
	| 'left'
	| 'letterSpacing'
	| 'lineHeight'
	| 'listStyleImage'
	| 'listStylePosition'
	| 'listStyleType'
	| 'listStyle'
	| 'marginRight'
	| 'marginLeft'
	| 'marginTop'
	| 'marginBottom'
	| 'margin'
	| 'maxHeight'
	| 'maxWidth'
	| 'minHeight'
	| 'minWidth'
	| 'orphans'
	| 'outlineColor'
	| 'outlineStyle'
	| 'outlineWidth'
	| 'outline'
	| 'overflow'
	| 'paddingTop'
	| 'paddingRight'
	| 'paddingBottom'
	| 'paddingLeft'
	| 'padding'
	| 'pageBreakAfter'
	| 'pageBreakBefore'
	| 'pageBreakInside'
	| 'pauseAfter'
	| 'pauseBefore'
	| 'pause'
	| 'pitchRange'
	| 'pitch'
	| 'playDuring'
	| 'position'
	| 'quotes'
	| 'richness'
	| 'right'
	| 'speakHeader'
	| 'speakNumeral'
	| 'speakPunctuation'
	| 'speak'
	| 'speechRate'
	| 'stress'
	| 'tableLayout'
	| 'textAlign'
	| 'textDecoration'
	| 'textIndent'
	| 'textTransform'
	| 'top'
	| 'unicodeBidi'
	| 'verticalAlign'
	| 'visibility'
	| 'voiceFamily'
	| 'volume'
	| 'whiteSpace'
	| 'widows'
	| 'width'
	| 'wordSpacing'
	| 'zIndex'
	| string