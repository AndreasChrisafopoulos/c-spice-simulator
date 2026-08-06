#include "parse.h"
#include "utils.h"
#include "stamp.h"

/** Performs the first pass of the netlist parsing.
 * This pass does not store elements; it only gathers information needed
 * for memory allocation and configuration.
 *
 * Specifically, it:
 * - Counts circuit elements and nodes
 * - Counts voltage sources and inductors (for MNA)
 * - Detects transient analysis (.tran) and transient sources
 * - Parses .options flags (sparse, spd, iter, custom, itol, method)
 */
void first_pass(
    const char *filename,   // netlist file name
    int *element_count,     
    int *node_count, 
    int *VsrcAndIndNumber, 
    int *use_cholesky, 
    int *use_custom,
    int *use_iterative,
    double *itol,
    int *use_sparse,
    int *has_trans,
    TRAN_Analysis *tran,
    int *tran_src_number 
) 
{
    int m2=0;
    int tran_src_cntr = 0;
    
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("open");
        exit(1);
    }

    NodeEntry *node_table = NULL;   // hash table head
    NodeEntry *n;
    int element_counter = 0, node_counter = 0;
    char line[MAX_LINE_LEN];

    tran->method = 1; // 1 = tr, default
    while (fgets(line, sizeof(line), fp)) {
        // Skip blank lines and comments
        char *p = line;
        while (isspace(*p)) p++;

        str_to_lower(p);
        if (strncasecmp(p, ".end", 4) == 0) break;

        // counts trans srcs
        if (strstr(line, "exp") || strstr(line, "sin") || strstr(line,"pwl") || strstr(line,"pulse"))
            tran_src_cntr++;

        // detect transient analysis 
        if (strncasecmp(p, ".tran", 5) == 0) {
            char cmd[8];
            double h, tstop;

            int parsed = sscanf(p, "%s %lf %lf", cmd, &h, &tstop);
            if(parsed != 3){
                printf("ERROR: Bad .TRAN format: %s\n", p);
                continue;
            }

            tran->time_step  = h;
            tran->final_time = tstop;
            tran->plot_count = 0;

            *has_trans = 1;
            continue;
        }

        if (strncasecmp(p, ".options", 8) == 0) {

            char *opt = p + 8;   // Part of the line after ".options"

            if (strcasestr(opt, "sparse"))
                *use_sparse = 1;

            // SPD  enables Cholesky or CG
            if (strcasestr(opt, "spd"))
                *use_cholesky = 1;

            // CUSTOM  your custom mode
            if (strcasestr(opt, "custom"))
                *use_custom = 1;

            // ITER  enables iterative methods
            if (strcasestr(opt, "iter"))
                *use_iterative = 1;

            if (strcasestr(opt, "method=tr"))
                tran->method = 1;
            if (strcasestr(opt, "method=be"))
                tran->method = 0;

            // ITOL=<value>
            char *itol_ptr = strcasestr(opt, "itol=");
            if (itol_ptr != NULL) {
                itol_ptr += 5;  // skip "itol="
                *itol = atof(itol_ptr);
            }
        }

        
        if (*p == '\0' || *p == '*' || *p == '.') continue;

        // Element type (R, C, L, V, I, D, M, Q)
        char type = toupper(*p);
        if (strchr("RCLVIDMQ", type) == NULL) continue;

        if (strchr("VL", type) !=NULL)
            m2++;
            
        element_counter++;

        // Tokenize
        char *tokens[10];
        int ntok = 0;
        char *tok = strtok(p, " \t\n");
        while (tok && ntok < 10) {
            str_to_lower(tok);
            tokens[ntok++] = tok;
            tok = strtok(NULL, " \t\n");
        }

        if (ntok < 3) continue;

        // Determine number of ports
        int ports = 0;
        if (type == 'R' || type == 'C' || type == 'L' ||
            type == 'V' || type == 'I' || type == 'D')
            ports = 2;
        else if (type == 'Q')
            ports = 3;
        else if (type == 'M')
            ports = 4;

        // Insert node names into hash
        for (int i = 1; i <= ports && i < ntok; i++) {
            HASH_FIND_STR(node_table, tokens[i], n);
            if (!n) {
                n = malloc(sizeof(NodeEntry));
                strncpy(n->name, tokens[i], sizeof(n->name));
                HASH_ADD_STR(node_table, name, n);
                node_counter++;
            }
        }
    }

    fclose(fp);
    *element_count = element_counter;
    *VsrcAndIndNumber = m2;
    *node_count = node_counter;
    *tran_src_number = tran_src_cntr;

    NodeEntry *curr, *tmp;
    HASH_ITER(hh, node_table, curr, tmp) {
        HASH_DEL(node_table, curr);
        free(curr);
    }
}

// Returns the enum type of an element.
elemenT getType(char p) 
{
    p = toupper(p); // ignore case
    switch (p) {
        case 'R': return RES;   // Resistor
        case 'C': return CAP;   // Capacitor
        case 'L': return IND;   // Inductor
        case 'V': return VSRC;  // Voltage Source
        case 'I': return ISRC;  // Current Source
        case 'D': return DIODE; // Diode
        case 'M': return MOS;   // MOS transistor
        case 'Q': return BJT;   // BJT transistor
        default:
            //fprintf(stderr, "Unknown element type: %c\n", p);
            return ZERO;
            //exit(EXIT_FAILURE);
    }
}

/** 
 * Creates a new node if it does not exist and returns its ID.
 * If the node already exists, its ID is returned.
 * This function also works as a find operation.
 * Assumes that netlist nodes used in DC, sweep, etc. already exist.
 */
int insert_node(NodeHashtable **table, const char *name, int table_size, int *node_counter) 
{
    unsigned int idx = hash_string(name, table_size);
    NodeHashtable *curr = table[idx];

    // Check if node already exists (case-insensitive)
    while (curr) {
        if (strcasecmp(curr->name, name) == 0)
            return curr->id;
        curr = curr->next;
    }

    // Create new node
    NodeHashtable *new_node = malloc(sizeof(NodeHashtable));
    strcpy(new_node->name, name);
    if (strcasecmp(name, "0") == 0) {
        new_node->id = 0;
    } else 
        new_node->id = (*node_counter)++;

    new_node->next = table[idx];

    // Insert into the linked list of this hash bucket
    table[idx] = new_node;

    return new_node->id;
}

/**
 * Maps an element name to its index in the elements array.
 *
 * Inserts a (name → index) entry into the ElementID hash table,
 * allowing constant-time lookup of elements stored in the main
 * elements array without performing O(n) searches.
 */
void insert_elementID(int index, const char *name, int table_size, ElementID **table) 
{
    unsigned int idx = hash_string(name, table_size);
    ElementID *curr = table[idx];

    // Traverse bucket (no duplicate check)
    while (curr) {
        curr = curr->next;
    }

    // Create new ElementID
    ElementID *new_element = malloc(sizeof(ElementID));
    strcpy(new_element->name, name);
    new_element->index = index;
    new_element->next = table[idx];

    // Insert into the linked list of this hash bucket
    table[idx] = new_element;
}


// Parses transient source specifications (EXP, SIN, PULSE, PWL) from a netlist line
// and fills the corresponding fields of the source structure.
// Returns 1 on success, 0 if no transient specification is found.
int parse_transient_spec(const char *line, TwoPortsSource *src)
{
    src->tr_type = TR_NONE;
    char *line_lower = line;
    char *start;

    str_to_lower(line_lower);
    //  EXP 
    start = strstr(line_lower, "exp");
    if(start){
        src->tr_type = TR_EXP;

        char *params = strchr(start, '(');   
        if(!params) return 0;

        sscanf(params,"(%lf %lf %lf %lf %lf %lf)",
               &src->tr.exp.i1,&src->tr.exp.i2,&src->tr.exp.td1,
               &src->tr.exp.tc1,&src->tr.exp.td2,&src->tr.exp.tc2);
        return 1;
    }

    //  SIN 
    start = strstr(line_lower, "sin");
    if(start){
        src->tr_type = TR_SIN;

        char *params = strchr(start,'(');
        if(!params) return 0;

        sscanf(params,"(%lf %lf %lf %lf %lf %lf)",
               &src->tr.sin.i1,&src->tr.sin.ia,&src->tr.sin.freq,
               &src->tr.sin.td,&src->tr.sin.df,&src->tr.sin.phase);
        return 1;
    }

    //  PULSE 
    start = strstr(line_lower,"pulse");
    if(start){
        src->tr_type = TR_PULSE;

        char *params = strchr(start,'(');
        if(!params) return 0;

        sscanf(params,"(%lf, %lf, %lf, %lf, %lf, %lf, %lf)",
               &src->tr.pulse.i1,&src->tr.pulse.i2,&src->tr.pulse.td,
               &src->tr.pulse.tr,&src->tr.pulse.tf,
               &src->tr.pulse.pw,&src->tr.pulse.per);
        return 1;
    }

    //  PWL 
    start = strstr(line_lower,"pwl");
    if(start){
        src->tr_type = TR_PWL;
        src->tr.pwl.points = 0;

        char *params = strchr(start,'(');
        while(params && src->tr.pwl.points < 20){
            double t,i;
            sscanf(params,"(%le %le)", &t,&i);

            int k = src->tr.pwl.points++;
            src->tr.pwl.t[k] = t;
            src->tr.pwl.i[k] = i;

            params = strchr(params+1,'('); 
        }
        return 1;
    }

    return 0;
}



/**
 * Parses a netlist element line, creates the corresponding element,
 * and inserts it into the circuit data structures.
 *
 * This function:
 * - Parses a single netlist line describing a circuit element
 * - Creates and initializes the appropriate element structure
 * - Inserts the element into the main elements or MNA-related element arrays
 * - Registers the element name in the ElementID hash table for O(1) lookup
 * - Creates or retrieves node IDs using a node name → ID hash table
 * - Stores node data in a sequential node array indexed by node ID
 * - Performs matrix stamping for DC and transient analysis
 * - Detects and registers transient-dependent sources when applicable
 *
 * Supported elements include:
 * R, C, L, V, I, D, M, Q (SPICE-like syntax).
 *
 * note: This function assumes that all tables and matrices have been
 *       pre-allocated by the caller.
 */

void insert_element(
    Element *elements,
    Element *m2_elements,
    int *m2_index, 
    int *index, 
    const char *line, 
    NodeHashtable **nodes_hashtable, 
    int nodeTable_size, 
    int elemntTable_size, 
    int *node_counter, 
    ElementID **elementIDs, 
    gsl_matrix *A, gsl_vector *b, 
    int *m2_index_forStamp, 
    Node *nodes,
    int use_sparse,
    cs *A_triplet,
    gsl_matrix *C,
    cs *C_triplet,
    int do_trans,
    int *m2_index_forStampC, TranSrc *tran_sources, int *tran_src_cntr
) 
{
    char name[MAX_NAME_LEN], n1[MAX_NAME_LEN], n2[MAX_NAME_LEN], n3[MAX_NAME_LEN], n4[MAX_NAME_LEN], model[MAX_NAME_LEN], area[MAX_NAME_LEN], l[MAX_NAME_LEN], w[MAX_NAME_LEN];
    double value;
    char *p = (char *)line;
    char t = toupper(*p);
    elemenT type = getType(*p);
    int elementsIDS_size = elemntTable_size;
    int id1, id2, id3, id4;
    
    // Dispatch element parsing, registration, and stamping based on element type
    if (strchr("RCL", t)) {
        sscanf(p, "%s %s %s %lf", name, n1, n2, &value);
        //TwoPortsElement *res;
        if(type == IND) {
            m2_elements[*m2_index].type = type;
            m2_elements[*m2_index].data = malloc(sizeof(TwoPortsElement));
            TwoPortsElement *res = (TwoPortsElement *)m2_elements[*m2_index].data;
            insert_elementID(*m2_index, name, elementsIDS_size, elementIDs);

            strcpy(res->name, name);
            res->value = value;
            res->ports[0] = insert_node(nodes_hashtable, n1, nodeTable_size, node_counter); //returns node ID, if its new id, stores it in node hashtable
            res->ports[1] = insert_node(nodes_hashtable, n2, nodeTable_size, node_counter);
            
            id1=res->ports[0];
            id2=res->ports[1];
            if (nodes[id1].id != id1) { //checks if the returned node already exists in node table
                nodes[id1].id = id1;    // stores the new node in node table
                strcpy(nodes[id1].name, n1); 
            }
            if (nodes[id2].id != id2) {
                nodes[id2].id = id2;
                strcpy(nodes[id2].name, n2); 
            }
            stamp(t, res->ports[0], res->ports[1], A, b, m2_index_forStamp, value, nodeTable_size/2, use_sparse, A_triplet);
            if (do_trans)
                stamp_C(t, res->ports[0], res->ports[1], C, m2_index_forStamp, value, nodeTable_size/2, use_sparse, C_triplet);
            (*m2_index)++;
        } else {
            elements[*index].type = type;
            elements[*index].data = malloc(sizeof(TwoPortsElement));
            TwoPortsElement *res = (TwoPortsElement *)elements[*index].data;
            insert_elementID(*index, name, elementsIDS_size, elementIDs);

            strcpy(res->name, name);
            res->value = value;
            res->ports[0] = insert_node(nodes_hashtable, n1, nodeTable_size, node_counter);
            res->ports[1] = insert_node(nodes_hashtable, n2, nodeTable_size, node_counter);
            
            id1=res->ports[0];
            id2=res->ports[1];
            if (nodes[id1].id != id1) {
                nodes[id1].id = id1;
                strcpy(nodes[id1].name, n1); 
            }
            if (nodes[id2].id != id2) {
                nodes[id2].id = id2;
                strcpy(nodes[id2].name, n2); 
            }
            stamp(t, res->ports[0], res->ports[1], A, b, m2_index_forStamp, value, nodeTable_size/2, use_sparse, A_triplet);
            if (do_trans)
                stamp_C(t, res->ports[0], res->ports[1], C, m2_index_forStamp, value, nodeTable_size/2, use_sparse, C_triplet);
            (*index)++;
        }
        
        
    } else if (strchr("VI", t)) {
        sscanf(p, "%s %s %s %lf", name, n1, n2, &value);
        if(type == VSRC) {
            m2_elements[*m2_index].type = type;
            m2_elements[*m2_index].data = malloc(sizeof(TwoPortsSource));
            TwoPortsSource *src = (TwoPortsSource *)m2_elements[*m2_index].data;
            insert_elementID(*m2_index, name, elementsIDS_size, elementIDs);

             strcpy(src->name, name);
            src->value = value;
            src->ports[0] = insert_node(nodes_hashtable, n1, nodeTable_size, node_counter);
            src->ports[1] = insert_node(nodes_hashtable, n2, nodeTable_size, node_counter);
            

            id1=src->ports[0];
            id2=src->ports[1];
            if (nodes[id1].id != id1) {
                nodes[id1].id = id1;
                strcpy(nodes[id1].name, n1); 
            }
            if (nodes[id2].id != id2) {
                nodes[id2].id = id2;
                strcpy(nodes[id2].name, n2); 
            }
            
            int exists = parse_transient_spec(line, src);
            if (exists) {
                tran_sources[*tran_src_cntr].index = *m2_index;
                tran_sources[*tran_src_cntr].type = 1;
                (*tran_src_cntr)++;
            }
            stamp(t, src->ports[0], src->ports[1], A, b, m2_index_forStamp, value, nodeTable_size/2, use_sparse, A_triplet);
            (*m2_index)++;
        } else {
            elements[*index].type = type;
            elements[*index].data = malloc(sizeof(TwoPortsSource));
            TwoPortsSource *src = (TwoPortsSource *)elements[*index].data;
            insert_elementID(*index, name, elementsIDS_size, elementIDs);

             strcpy(src->name, name);
            src->value = value;
            src->ports[0] = insert_node(nodes_hashtable, n1, nodeTable_size, node_counter);
            src->ports[1] = insert_node(nodes_hashtable, n2, nodeTable_size, node_counter);
            

            id1=src->ports[0];
            id2=src->ports[1];
            if (nodes[id1].id != id1) {
                nodes[id1].id = id1;
                strcpy(nodes[id1].name, n1); 
            }
            if (nodes[id2].id != id2) {
                nodes[id2].id = id2;
                strcpy(nodes[id2].name, n2); 
            }
            
            int exists = parse_transient_spec(line, src);
            if (exists) {
                tran_sources[*tran_src_cntr].index = *index; //    [*tran_src_cntr][0] = *index;
                tran_sources[*tran_src_cntr].type = 0;    //[1] = 0;
                (*tran_src_cntr)++;
            }
            stamp(t, src->ports[0], src->ports[1], A, b, m2_index_forStamp, value, nodeTable_size/2, use_sparse, A_triplet);
            (*index)++;
        }
       
    } else if (t == 'D') {
        int count = sscanf(p, "%s %s %s %s %s", name, n1, n2, model, area);
        elements[*index].type = DIODE;
        elements[*index].data = malloc(sizeof(Diode));
        Diode *d = (Diode *)elements[*index].data;
        strcpy(d->name, name);
        strcpy(d->model_name, model);
        d->ports[0] = insert_node(nodes_hashtable, n1, nodeTable_size, node_counter);
        d->ports[1] = insert_node(nodes_hashtable, n2, nodeTable_size, node_counter);
        d->area = 1.0;
        if (count == 5 && strncmp(area, "area=", 5) == 0)
            d->area = atof(area + 5);
        insert_elementID(*index, name, elementsIDS_size, elementIDs);

        id1=d->ports[0];
        id2=d->ports[1];
        if (nodes[id1].id != id1) {
            nodes[id1].id = id1;
            strcpy(nodes[id1].name, n1); 
        }
        if (nodes[id2].id != id2) {
            nodes[id2].id = id2;
            strcpy(nodes[id2].name, n2); 
        }

        (*index)++;
    } else if (t == 'M') {
        sscanf(p, "%s %s %s %s %s %s %s %s", name, n1, n2, n3, n4, model, l, w);
        elements[*index].type = MOS;
        elements[*index].data = malloc(sizeof(Mos));
        Mos *m = (Mos *)elements[*index].data;
        strcpy(m->name, name);
        strcpy(m->model_name, model);
        m->l = atof(l + 2);
        m->w = atof(w + 2);
        m->ports[0] = insert_node(nodes_hashtable, n1, nodeTable_size, node_counter);
        m->ports[1] = insert_node(nodes_hashtable, n2, nodeTable_size, node_counter);
        m->ports[2] = insert_node(nodes_hashtable, n3, nodeTable_size, node_counter);
        m->ports[3] = insert_node(nodes_hashtable, n4, nodeTable_size, node_counter);
        insert_elementID(*index, name, elementsIDS_size, elementIDs);

        id1=m->ports[0];
        id2=m->ports[1];
        id3=m->ports[2];
        id4=m->ports[3];
        if (nodes[id1].id != id1) {
            nodes[id1].id = id1;
            strcpy(nodes[id1].name, n1); 
        }
        if (nodes[id2].id != id2) {
            nodes[id2].id = id2;
            strcpy(nodes[id2].name, n2); 
        }
        if (nodes[id3].id != id3) {
            nodes[id3].id = id3;
            strcpy(nodes[id3].name, n3); 
        }
        if (nodes[id4].id != id4) {
            nodes[id4].id = id4;
            strcpy(nodes[id4].name, n4); 
        }
        (*index)++;
    } else if (t == 'Q') {
        int count = sscanf(p, "%s %s %s %s %s %s", name, n1, n2, n3, model, area);
        elements[*index].type = BJT;
        elements[*index].data = malloc(sizeof(Bjt));
        Bjt *q = (Bjt *)elements[*index].data;
        strcpy(q->name, name);
        strcpy(q->model_name, model);
        q->ports[0] = insert_node(nodes_hashtable, n1, nodeTable_size, node_counter);
        q->ports[1] = insert_node(nodes_hashtable, n2, nodeTable_size, node_counter);
        q->ports[2] = insert_node(nodes_hashtable, n3, nodeTable_size, node_counter);
        q->area = 1.0;
        if (count == 6 && strncmp(area, "area=", 5) == 0)
            q->area = atof(area + 5);

        id1=q->ports[0];
        id2=q->ports[1];
        id3=q->ports[2];
        if (nodes[id1].id != id1) {
            nodes[id1].id = id1;
            strcpy(nodes[id1].name, n1); 
        }
        if (nodes[id2].id != id2) {
            nodes[id2].id = id2;
            strcpy(nodes[id2].name, n2); 
        }
        if (nodes[id3].id != id3) {
            nodes[id3].id = id3;
            strcpy(nodes[id3].name, n3); 
        }
        insert_elementID(*index, name, elementsIDS_size, elementIDs);
        (*index)++;
    }
}

/**
 * Parses a SPICE-like netlist file and builds the circuit data structures.
 *
 * This function:
 * - Reads the netlist line by line
 * - Handles analysis directives (.DC, .TRAN)
 * - Registers plot/print requests for DC sweep and transient analysis
 * - Delegates element parsing and stamping to insert_element()
 * - Builds node and element lookup tables during parsing
 *
 * note: Assumes that memory for elements, nodes, matrices, and sweep
 *       structures has already been allocated by the caller.
 */
void parse_netlist(
    const char *filename, 
    Element *elements, 
    Element *m2_elements, 
    NodeHashtable **nodes_hashtable, 
    int nodeTable_count, 
    int elementTable_size, 
    int *node_counter, 
    ElementID **elementIDs, 
    gsl_matrix *A, 
    gsl_vector *b, 
    int *m2_index_forStamp, 
    Node *nodes, 
    DCSweep *dc_sweeps, 
    int *dc_sweep_count,
    int use_sparse,
    cs *A_triplet,
    gsl_matrix *C,
    cs *C_triplet,
    int do_trans,
    TRAN_Analysis *tran, TranSrc *tran_sources
) 
{   
    int m2_index_forStampC=0;
    int in_dc_section=0;
    int in_tran_section=0;
    int tran_src_cntr = 0;


    int i = 0, m2_index=0;
    *node_counter = 1;
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("File open failed");
        exit(EXIT_FAILURE);
    }

    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        char *p = line;
        while (isspace(*p)) p++;

        //str_to_lower(p);
        if (strncasecmp(p, ".end", 4) == 0) break;

        if (strncasecmp(p, ".dc", 3) == 0) {
            in_dc_section = 1;
            in_tran_section = 0;
            char cmd[8], var[64];
            double start, end, step;

            int parsed = sscanf(p, "%s %s %lf %lf %lf", cmd, var, &start, &end, &step);

            if (parsed != 5) {
                printf("ERROR: Bad .DC format: %s\n", p);
                continue;
            }

        
            strcpy(dc_sweeps[*dc_sweep_count].source_name, var);
            dc_sweeps[*dc_sweep_count].start = start; 
            dc_sweeps[*dc_sweep_count].end = end;
            dc_sweeps[*dc_sweep_count].step = step;
            dc_sweeps[*dc_sweep_count].plot_count = 0;
            (*dc_sweep_count)++;
        }

       
        //  TRAN detect 
        if(strncasecmp(p,".tran",5)==0){
            in_dc_section = 0;
            in_tran_section = 1;
            continue;
        }

        // PLOT / PRINT 
        // PLOT / PRINT 
    if(!strncasecmp(p,".plot",5) || !strncasecmp(p,".print",6)){

        char cmd[8], rest[256];
        sscanf(p, "%s %[^\n]", cmd, rest);
        char *token = strtok(rest, " ");

        // Plot DC 
        if(in_dc_section){
            while(token){
                if(dc_sweeps[*dc_sweep_count-1].plot_count >= MAX_SWEEP_PLOTS)
                    break;

                if( (token[0]=='V' || token[0]=='v') && token[1]=='(' ){
                    char node[64];
                    sscanf(token, "%*[Vv](%[^)])", node);

                    int id = insert_node(
                        nodes_hashtable,
                        node,
                        nodeTable_count,
                        node_counter
                    );

                    dc_sweeps[*dc_sweep_count-1].plot_nodes[
                        dc_sweeps[*dc_sweep_count-1].plot_count++
                    ] = id;
                }

                token = strtok(NULL, " ");
            }
            continue;
        }

        // Plot TRAN 
        if(in_tran_section){
            while(token){
                if(tran->plot_count >= MAX_SWEEP_PLOTS)
                    break;

                if( (token[0]=='V' || token[0]=='v') && token[1]=='(' ){
                    char node[64];
                    sscanf(token, "%*[Vv](%[^)])", node);

                    int id = insert_node(
                        nodes_hashtable,
                        node,
                        nodeTable_count,
                        node_counter
                    );

                    tran->plot_nodes[tran->plot_count++] = id;
                }

                token = strtok(NULL, " ");
            }
            continue;
        }
    }
            

        if (*p == '*' || *p == '\n' || *p == '.' || *p == '\0') // ignore comments / empty lines / directives
            continue;

        insert_element(elements, m2_elements, &m2_index, &i, p, nodes_hashtable, nodeTable_count, elementTable_size, node_counter, elementIDs, A, b, m2_index_forStamp, nodes, use_sparse, A_triplet, C, C_triplet, do_trans, &m2_index_forStampC, tran_sources, &tran_src_cntr);
    }
    fclose(fp);
}





