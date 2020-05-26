import nose
from nose.tools import assert_equals
import smcsmc

# This has real data for now, but I will change it over later once it is working.
def test_init():
    o = smcsmc.Output("/users/lunter/ccole/repos/ancient_migration/analyses/full_africa/S_Yoruba-1_10000/result.out")
