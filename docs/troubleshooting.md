# Troubleshooting
## Bytecode error
On some systems, the installed version of Java thinks there is an error in the Java bytecode. If this is the case, 
you will probably get an error somewhat like this:
```
Exception in thread "main" java.lang.VerifyError: Bad <init> method call from inside of a branch
Exception Details:
  Location:
    graxxia/Matrix.<init>(Ljava/lang/Iterable;)V @151: invokespecial
  Reason:
    Error exists in the bytecode
  Bytecode:
    0x0000000: b800 2b4d 05bd 002d 5903 2b12 06b8 0092
    0x0000010: 5359 0401 129e b800 9253 5910 ff12 02b8
    0x0000020: 00bb 2a5f ab00 0000 0000 0169 0000 000b
    0x0000030: 942a 6887 0000 0064 a152 4770 0000 0079
    0x0000040: b12c cb16 0000 008e c4ba 1d9e 0000 00a3
    0x0000050: dd06 2917 0000 00c2 ddc5 4bfe 0000 00d3
    0x0000060: e9fe b586 0000 00e8 0255 f295 0000 010b
    0x0000070: 2b48 1cef 0000 0122 496e 9609 0000 0143
    0x0000080: 74b6 5b2d 0000 0154 5f5a 5903 3212 bdb8
    0x0000090: 0048 c000 bd5f 57b7 00bf a700 fd5f 5a59
    0x00000a0: 0332 12c1 b800 48c0 00c1 5f57 b700 c3a7
    0x00000b0: 00e8 5f5a 5903 3212 c5b8 0048 c000 c55f
    0x00000c0: 57b7 00c7 a700 d35f 5a59 0332 b800 cb5f
    0x00000d0: 5904 32b8 00cb 5f05 12b4 b800 cfc0 00b4
    0x00000e0: b700 d1a7 00b4 5f5a 0312 d2b8 00cf c000
    0x00000f0: d2b7 00d4 a700 a35f 5a59 0332 1206 b800
    0x0000100: 48c0 0006 5f57 b700 d6a7 008e 5f5a 5903
    0x0000110: 32b8 00cb 5f59 0432 b800 cb5f 5905 3212
    0x0000120: 9eb8 0048 c000 9e5f 57b7 00d8 a700 6b5f
    0x0000130: 5a59 0332 b800 cb5f 5904 32b8 00cb 5f57
    0x0000140: b700 daa7 0054 5f5a 5903 3212 06b8 0048
    0x0000150: c000 065f 5904 3212 9eb8 0048 c000 9e5f
    0x0000160: 57b7 00dc a700 335f 5a03 128b b800 cfc0
    0x0000170: 008b b700 dea7 0022 5f5a 5903 3212 38b8
    0x0000180: 0048 c000 385f 57b7 00e0 a700 0dbb 00e2
    0x0000190: 5912 e4b7 00e7 bf57 b1                 
  Stackmap Table:
    full_frame(@136,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@157,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@178,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@199,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@230,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@247,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@268,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@303,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@326,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@359,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@376,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@397,{UninitializedThis,Object[#6],Object[#160]},{Object[#233],UninitializedThis})
    full_frame(@407,{Object[#2],Object[#6],Object[#160]},{Object[#233]})
```

This isn't a significant error, so you can make Java stop reporting it by editing your `config.groovy` and
setting `JAVA_OPTS="-noverify"`.

## Out of Memory
If you notice OutOfMemory errors in the cpipe output log, especially when you are running a large number of samples, you
might need to increase the amount of memory allocated to Java. You can do this by setting the `MAX_JAVA_MEM` variable
like this:
```bash
MAX_JAVA_MEM=8g cpipe run
```

