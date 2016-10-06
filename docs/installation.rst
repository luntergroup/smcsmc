.. _sec-installation:

============
Installation
============

``smcsmc`` is written in C++.

**************
Stable Release
**************

.. todo::
    You can download the latest stable release packaged for a variety of different platform from


*******************************
Development Version From GitHub
*******************************

You can also install ``smcsmc`` directly from the git repository. Here, you will need ``autoconf``, check whether this is already installed by running:

.. code-block:: bash

    $ which autoconf

On Debian/Ubuntu based systems:

.. code-block:: bash

    $ apt-get install build-essential autoconf autoconf-archive libcppunit-dev


On Mac OS:

.. code-block:: bash

    $ port install automake autoconf autoconf-archive cppunit


Afterwards you can clone the code from the github repository,

.. code-block:: bash

    $ git clone git@github.com:luntergroup/smcsmc.git
    $ cd smcsmc

and build the binary using

.. code-block:: bash

    $ ./bootstrap
    $ make

or install with ``make install``.
