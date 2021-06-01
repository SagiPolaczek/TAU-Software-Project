import capi_demo2 as cd2
# To use this run in the terminal: python setup.py build_ext --inplace
# And on Nova: python3.8.5 setup.py build_ext --inplace
cd2.f1()
print(cd2.f3([9,9]))
print(cd2.f3((9,100)))
cd2.f4(5, "second", "third")  # v4 is not read and keeps his value
# cd2.f4(5,"A", "B", [33])  # TypeError: What are you doing?
cd2.raise_hell(4)
# See how in the TypeError it says f3()
print(cd2.f3(9,9))