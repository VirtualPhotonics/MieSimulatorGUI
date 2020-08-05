When contributing to the Mie Simulator GUI please follow the coding conventions detailed below. Some of the initial development did not adhere to these guidelines but going forward, this will be the code style and structure.
## Naming Conventions
### Definitions:
* __Camel Case:__ Naming convention where the first letter of each word is capitalized with the exception of the first letter. Example: __thisIsCamelCased__
* __Pascal Case:__ Naming convention where the first letter of each word is capitalized including the first letter. Example: __ThisIsPascalCased__
* __Underscore Uppercase Delimited:__ Naming convention that places underscores between the words and the first letter of each word is capitalized. Example: __This_Is_Underscore_Delimited__
* __Underscore Lowercase Delimited:__ Naming convention that places underscores between the words and the first letter of each word is not capitalized. Example: __this_is_underscore_lowercase_delimited__<br />

These are the only naming conventions that should be used in the Mie Simulator GUI and apply to different aspects of the code defined below.
* __Classes__, __properties__, and __functions__, should be Pascal cased. Example: __MyClass__
* __Private variables__  (except UI variables) should be prefixed with an underscore and camel cased. Example: ___privateVariables__
* __Member variables__, __parameters__ , __UI variables__, and __local variables__</strong> should be camel cased. Example: __memberVariable__
* __UI Methods__ should be Underscore Lowercase Delimited. Example: __on_object_name_action__
* __Test methods__ should be underscore delimited. Example: __Test_Methods__



## Indenting
Code should be indented with 4 space characters, if a tab character is used, make sure your source-code editor replaces tabs with 4 spaces.
## Braces
An opening brace { should appear on the line after the start of a statement block and the code inside the brace should be indented 4 spaces. The closing brace } should be inline with the opening brace on it's own line. 

Braces may may be avoided when there is only a single line of code.  

## Spacing
Spacing improves the readability of the source-code, follow these guidlines for spacing:<br /> 

 * Use a space after the comma between arguments in a function
 * Use a space after statements like __if__, __while__, and __for
 * Use a space before and after operators with the exception of ++ and --

Example:

```c++
#include <math.h>

class Complex
{
public:
    Complex(double re, double im)
        : _re(re), _im(im)
    {}
    double modulus() const
    {
        return sqrt(_re * _re + _im * _im);
    }
private:
    double _re;
    double _im;
};

void bar(int i)
{
    static int counter = 0;
    counter += i;
}

namespace Foo
{
namespace Bar
{
void foo(int a, int b)
{
    for (int i = 0; i < a; i++)
    {
        if (i < b)
            bar(i);
        else
        {
            bar(i);
            bar(b);
        }
    }
}
} // namespace Bar
} // namespace Foo
```
## Comments
The style for the code comments is two slashes // even for multi-line comments. Avoid using /* comments */

