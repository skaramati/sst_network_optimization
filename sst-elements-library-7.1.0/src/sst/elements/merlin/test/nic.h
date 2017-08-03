// -*- mode: c++ -*-

// Copyright 2009-2017 Sandia Corporation. Under the terms
// of Contract DE-NA0003525 with Sandia Corporation, the U.S.
// Government retains certain rights in this software.
// 
// Copyright (c) 2009-2017, Sandia Corporation
// All rights reserved.
// 
// Portions are copyright of other developers:
// See the file CONTRIBUTORS.TXT in the top level directory
// the distribution for more information.
//
// This file is part of the SST software package. For license
// information, see the LICENSE file in the top level directory of the
// distribution.


#ifndef COMPONENTS_MERLIN_TEST_NIC_H
#define COMPONENTS_MERLIN_TEST_NIC_H

#include <sst/core/component.h>
#include <sst/core/event.h>
#include <sst/core/elementinfo.h>
#include <sst/core/link.h>
#include <sst/core/timeConverter.h>
#include <sst/core/interfaces/simpleNetwork.h>


namespace SST {

namespace Merlin {


class nic : public Component {

public:

    SST_ELI_REGISTER_COMPONENT(
        nic,
        "merlin",
        "test_nic",
        SST_ELI_ELEMENT_VERSION(1,0,0),
        "Simple NIC to test base functionality.",
        COMPONENT_CATEGORY_NETWORK)
    
    SST_ELI_DOCUMENT_PARAMS(
        {"id",           "Network ID of endpoint."},
        {"num_peers",    "Total number of endpoints in network."},
        {"num_messages", "Total number of messages to send to each endpoint."},
        {"num_vns",      "Number of requested virtual networks."},
        {"link_bw",      "Bandwidth of the router link specified in either b/s or B/s (can include SI prefix)."},
        {"topology",     "Name of the topology subcomponent that should be loaded to control routing."},
        {"remap",        "Creates a logical to physical mapping shifted by remap amount.", "0"}
    )

    SST_ELI_DOCUMENT_PORTS(
        {"rtr",  "Port that hooks up to router.", { "merlin.RtrEvent", "merlin.credit_event" } }
    )


private:

    // SST::Interfaces::SimpleNetwork::nid_t id;
    int id;
    int net_id;
    int num_peers;
    int num_vns;
    int last_vn;

    int num_msg;
    int packets_sent;
    int packets_recd;
    int stalled_cycles;

    bool done;
    bool initialized;
    
    SST::Interfaces::SimpleNetwork* link_control;

    int last_target;
    
    int *next_seq;

    int remap;

    Output& output;
    
public:
    nic(ComponentId_t cid, Params& params);
    ~nic();

    void init(unsigned int phase);
    void setup(); 
    void finish();


private:
    bool clock_handler(Cycle_t cycle);

};

}
}

#endif // COMPONENTS_MERLIN_TEST_NIC_H
