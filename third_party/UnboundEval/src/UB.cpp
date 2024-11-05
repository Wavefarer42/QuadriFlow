
#define UUID_SYSTEM_GENERATOR
#include "uuid.h"

#include "UB.h"

uuids::uuid UB::generateUUID()
{
    return uuids::uuid_system_generator{}();
}
