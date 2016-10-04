const DEFAULT_INDENT = '  '
const KEYWORD_REGEXP = new RegExp
('^('+
  'abstract|boolean|break|byte|case|catch|char|class|const|continue|debugger|'+
  'default|delete|do|double|else|enum|export|extends|false|final|finally|'+
  'float|for|function|goto|if|implements|import|in|instanceof|int|interface|'+
  'long|native|new|null|package|private|protected|public|return|short|static|'+
  'super|switch|synchronized|this|throw|throws|transient|true|try|typeof|'+
  'undefined|var|void|volatile|while|with'+
')$')


var _filter
var _indent = DEFAULT_INDENT

var seen = []


function applyToSource(prototype, func)
{
  if(prototype.toSource === undefined)
    Object.defineProperty(prototype, 'toSource', {value: func})
}

function fill(keys, func, thisArg)
{
  thisArg = thisArg || keys

  var items = ''

  if(keys.length)
  {
    seen.push(thisArg)

    items = join(keys.map(func, thisArg))

    seen.pop()
  }

  return items
}

function getIndent(indent)
{
  switch(typeof indent)
  {
    case 'boolean':   return indent ? DEFAULT_INDENT : ''
    case 'number':    return Array(indent+1).join(' ')
    case 'string':    return indent
    case 'undefined': return DEFAULT_INDENT
  }

  if(indent === null) return ''

  throw SyntaxError('Invalid indent: '+indent)
}

function join(elements)
{
  if(elements.length < 2) return elements

  var offset = Array(seen.length).join(_indent)

  var newLine = (_indent&&'\n')+offset+_indent

  return newLine
       + elements.join(','+newLine)
       + (_indent&&'\n')+offset;
}

function legalKey(string)
{
  return /^[a-z_$][0-9a-z_$]*$/gi.test(string) && !KEYWORD_REGEXP.test(string)
}

function walk(object)
{
  object = _filter ? _filter(object) : object

  if (object === null) return 'null'
  if (typeof object === 'undefined') return 'undefined'

  var index = seen.indexOf(object);
  if (index >= 0) return '{$circularReference:'+index+'}'

  return object.toSource();
}


// Objects where toSource is equal to toString

[Boolean, Function, Number, RegExp].forEach(function(constructor)
{
  var prototype = constructor.prototype

  applyToSource(prototype, prototype.toString)
});


// Regular objects

applyToSource(Date.prototype, function()
{
  return 'new Date('+this.getTime()+')'
})

applyToSource(String.prototype, function()
{
  return JSON.stringify(this)
})

applyToSource(Math, function()
{
  return 'Math'
})


// Recursive objects (Array & Object)

applyToSource(Array.prototype, function()
{
  return '[' + fill(this, walk) + ']'
})

applyToSource(Object.prototype, function()
{
  var keys = Object.keys(this)

  function keyValue(key) {
    var s_key = legalKey(key) ? key : JSON.stringify(key)
    var value = walk(this[key])

    return s_key + ': ' + value
  }

  return '{' + fill(keys, keyValue, this) + '}'
})

//
// module.exports = function(object, filter, indent)
// {
//   _filter = filter
//   _indent = getIndent(indent)
//
//   var result = walk(object)
//
//   _filter = undefined
//   _indent = DEFAULT_INDENT
//
//   return result
// }
