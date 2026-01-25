# Contribution Workflow

All contributions to the **MieSimulatorGUI** must follow the workflow outlined below to ensure changes are traceable, reviewable, and consistent.

## 1. Create an Issue First
Before starting any work:
* **Create a GitHub issue** describing the bug, enhancement, or feature request.
* Clearly explain the motivation, scope, and expected behavior.
* Wait for confirmation or discussion if the change is significant.
* *Every code change must be associated with an issue.*

## 2. Fork the Repository
* Fork the main repository into your own GitHub account.
* All development work must be done in your fork, not directly in the main repository.

## 3. Create a Feature Branch
Create a new branch from the default branch (e.g., `master`). The branch name must be derived from the issue and prepended with `feature/`.

**Branch naming format:**
`feature/<issue-number>-short-description`

**Example:**
`feature/42-add-particle-size-slider`

## 4. Implement the Change
* Follow all coding conventions defined in this document.
* Keep changes focused on the scope of the issue.
* Ensure the code builds and tests pass before submitting.

## 5. Open a Pull Request
Open a Pull Request (PR) from your feature branch to the main repository. Reference the related issue in the PR description (e.g., `Closes #42`).

**Clearly describe:**
* What was changed
* Why the change was made
* Any relevant design or implementation details

## 6. Code Review and Approval
* All pull requests require review before merging.
* Be prepared to address feedback or requested changes.
* Once approved, the pull request will be merged by a repository maintainer.

> [!IMPORTANT]
> Direct commits to the main repository are not permitted.

# Coding Convention
When contributing to the **MieSimulatorGUI**, please follow the coding conventions detailed below. These guidelines represent the established style observed in the core simulation engine (e.g., `Calculate` and `MieSimulation` classes).

## 1. Naming Conventions

### Definitions:
* **Camel Case:** Naming convention where the first letter of each word is capitalized with the exception of the first letter. Example: `thisIsCamelCase`
* **Pascal Case:** Naming convention where the first letter of every word is capitalized. Example: `ThisIsPascalCase`
* **M-Prefix Pascal Case:** Member variables prefixed with a lowercase 'm' followed by the variable name in Pascal Case. Example: `mMemberVariable`

### Application:
* **Source and Header Files:** Must be **lowercase**.
  * *Example:* `miesimulation.cpp`, `miesimulation.h`, .
* **Classes and Functions:** Must be **Pascal Case**. 
  * *Example:* `class Calculate`, `void DoSimulation()`, `double CalculateG()`.
* **Member Variables:** Must use **M-Prefix Pascal Case**. 
  * *Example:* `mWavel`, `mQSca`, `mMinTheta`.
* **Local Variables and Parameters:** Must be **Camel Case**. 
  * *Example:* `stepTheta`, `para`, `curMus`.
* **UI Methods (Qt Slots):** Should follow the standard Qt underscore convention. 
  * *Example:* `on_button_clicked`.
* **UI Elements :** Should follow the underscore with Pascal Case. 
  * *Example:* `QRadioButton *radioButton_NumDen`.
* **Constants and Macros:** Must be **Uppercase Underscore Delimited**. 
  * *Example:* `CALCULATE_H`, `M_PI`.

## 2. Indenting
Code should be indented with 4 space characters, if a tab character is used, make sure your source-code editor replaces tabs with 4 spaces.

## 3. Braces
An opening brace { should appear on the line after the start of a statement block and the code inside the brace should be indented 4 spaces. The closing brace } should be inline with the opening brace on it's own line.<br /> 

## 4. Spacing
Spacing improves the readability of the source-code, follow these guidlines for spacing:
* Use a **space after commas** between function arguments: `void Function(int a, int b)`.
* Use a **space after control statements** like `if`, `while`, and `for`: `for (unsigned int i = 0; ...)`.
* Use a **space before and after binary operators**: `x = a + b` or `refRelRe = para->scatRefRealArray[r] / para->medRefArray[r]`.
* **No space** between a function name and its opening parenthesis: `CalculateG(curS1, curS2, para)`.

## 5. Example
The following snippet demonstrates the conventions used in the project:

```c++
#ifndef SAMPLE_H
#define SAMPLE_H

class Simulation
{
public:
    Simulation(void);
    
    double mResultValue; // Member variable with 'm' prefix

    void RunProcess(Parameters *para)
    {
        double localStep = 0.1; // Local variable in Camel Case
        
        for (int i = 0; i < para->nTheta; i++)
        {
            if (i % 2 == 0)
                CalculateInternal(i);
            else
            {
                mResultValue += i * localStep;
                UpdateLog();
            }
        }
    }

private:
    void CalculateInternal(int value);
};

#endif

