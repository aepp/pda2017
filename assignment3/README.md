How to:
1. Execute ``make clean && make`` in root directory to link and compile
2. Run with ``mpiexec -np [number of processes] ./assignment [program parameters]``
3. Enjoy!

Notes:

Ich habe die Stellen mit ``-->`` markiert, an denen ich vorher dynamisch den Speicher allokiert habe. Dies hat für 4 und 8 Prozesse funktioniert, bei 1 und 2 gab's dann Fehlermeldungen. Ich weiß, dass man den auf diese Weise allokierten Speiche auch mit ``free()`` freigeben muss, jedoch hat ``free()`` ebenfalls eine Fehlermeldung ausgespuckt, dass der übergeben Zeiger ungültig sei.