function logo()
txt1 = """
                   |
 ___  ___  ___   __|
(  ,\\(  _)(  ,) (  )|
 ) _/ ) _) )  \\  )(|
(_)  (___)(_)\\_)(__)|
                    |
"""

txt2 = """
 ___  _  _  _  _ |
(   \\( \\/ )( \\( )|
) ) )\\  /  )  ( |
(___/(__/  (_)\\_)|
           |
           |
"""

txt3 = """
           |
           |
 __  __   |
(  )(  )  |
 _  __)(  )(__ |
(_)(___/ (____)|
"""

fn(txt, mcro) = mapreduce((x) -> mcro(x)*"\n", *,
            split(txt, '\n'))

# height = 8
# a = TextBox(fn(txt1, (x)->@red(x)), width=22, padding=(0, 0, 0, 0), height=height)
# b = TextBox(fn(txt2, (x)->@green(x)), width=19, padding=(0, 0, 0, 0), height=height)
# c = TextBox(fn(txt3, (x)->@magenta(x)), width=17, padding=(0, 0, 0, 0), height=height)
# println(Term.hstack(a, b, c; pad=0))

for (l1,l2,l3) in zip(split(txt1, "|\n"),
                    split(txt2, "|\n"),
                    split(txt3, "|\n"))
    print(@red(@bold l1), @green(@bold l2), @magenta(@bold l3))
    println()
end
version()
end

logo()
