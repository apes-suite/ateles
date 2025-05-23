# Contributing to Ateles

Thanks for your interest in the development.
When contributing code, please heed the
[Coding Guidelines](https://apes-suite.github.io/treelm//page/codingGuidelines.html).

The project is organized in submodules, due to the
use of various components in different tools of the
APES framework.
You find the actual source code in the
[ateles-source](https://github.com/apes-suite/ateles-source)
repository.
To ease the work with the submodules, there is a
`request` script in the bin submodule.
It requires the [github cli](https://cli.github.com/)
and takes care of creating appropiate pull-requests.
For the interaction with the repository you should
always put you your work onto a branch.

Use this `request` script instead of `git push`, when
you want to put something back to the github
repository.

Thus, in general the workflow looks akin to this:

```
cd atl
git checkout -b my-work-branch
# modify code
../bin/request
```

The first invocation of `request` on the branch will
create a new pull request, and also create one in the
parent repository.
Subsequent calls will push your commits and update the
parent repository.
