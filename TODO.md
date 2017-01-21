# TODO

* fix luneburg function for the sporadic values it does not work, such as
    * 2^21, 2^30, 2^31
    * 3^12
    * 5^16
    * 7^12
    * 11^15
    * 13^12

Le cas 2^31 a l'air d'être le plus simple à étudier. Ce serait plus simple avec
Julia qui fait les calculs à notre place. D'autant plus que Julia est basé sur
Flint.
