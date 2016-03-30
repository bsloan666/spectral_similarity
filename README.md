# spectral_similarity
Python implementation of provisional SSI

Once unpacked it should be possible to enter the top-level source directory and type:

python ./test_illuminant.py resources/emitters/json/generic_led_cool.json

or 

python ./test_illuminant.py resources/emitters/csv/LED_AastroLightRing_3320K.csv

This will print a number representing the SSI (0-100 scale) for the test spectrum provided as the program's first commandline argument

The test is currently hard-wired to use ISO 7589 Tungsten as the reference illuminant.



