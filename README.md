# aaalligator
AAA Rational Function Approximation with rust

An attempt at implementing a rust variant of the methods described in "AAA rational approximation on a continuum" 2023 (Toby Driscoll, Yuji Nakatsukasa, Lloyd N. Trefethen) 10.48550/arXiv.2305.03677 https://arxiv.org/abs/2305.03677

Forgive me am i fairly new to programming so i'm pretty sure this code is full of bugs and things i can't work out how to fix yet. It's also missing a fair few features mentioned in the paper like lawson iterations and domains other than [-1, 1] and there are a few logic issues in the code but it otherwise appears to mostly work! :)

I ended up having to make two variants of the code aaaxf (floating point functions) and aaaxc (complex functions) due to not knowing if there was a way to force rust to take either as an input like matlab and julia both often seem to have automatic type conversion.

Weirdly i seem to have issues running this program in debug (non --release mode) that seems to cause an overflow but only on windows and not in Linux but can be fixed by turning opt-level = 1 for the dev profile.

There are several usage examples with poles and zeros generated however i think the actual representation should be returned slightly differently i have more work to do to make it do the same as other packages this is just what i have so far :|

No i have not worked out how to make the returned function r callable yet too much rust lifetimes making things confusing :\
