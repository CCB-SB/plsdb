#ORIGINAL CODE TAKEN FROM MASTER THESIS OF MARKUS DILLMANN
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import pronto, json, re
import nltk
nltk.download('punkt')
nltk.download('wordnet')
from nltk.stem import WordNetLemmatizer
from nltk.tokenize import word_tokenize



class Disease_ontology:
    def __init__(self, location="doid.owl"):
        self.ont = pronto.Ontology(location)
        self.relation_dict = dict()
        self.build_graph_doid()

    def build_graph_doid(self):
        """
        parsed the owl into graph
        """
        for term in self.ont:            
            # step over deprecated data
            if "deprecated" in term.other:
                if term.other["deprecated"]:
                    continue
            # remove some words like cell or cells
            stopwords = []
            querywords = term.name.split()
            resultwords = [word for word in querywords if word.lower() not in stopwords]
            term.name = ' '.join(resultwords)
            term.name = term.name.lower()

            self.relation_dict[term.name] = {'id': term.id, "can_be": {}, "is_a": {},
                                             "hasRelatedSynonym": set(), "related_to": {}, "is_a_childs": set(),
                                             "hasExactSynonym": set(),
                                             'xref': set(), 'has_symptom': [], 'derives_from': []}

            for key, value in term.relations.items():
                querywords = value.name[0].split()
                resultwords = [word for word in querywords if word.lower() not in stopwords]
                value.name[0] = ' '.join(resultwords)
                value.name[0] = value.name[0].lower()

                relation = str(key).split("'")[1].split("'")[0]
                if "name" not in self.relation_dict[term.name][relation] and "id" not in self.relation_dict[term.name][
                    relation]:
                    self.relation_dict[term.name][relation]["name"] = value.name
                    self.relation_dict[term.name][relation]["id"] = value.id
                else:
                    self.relation_dict[term.name][relation]["name"].append(value.name[0])
                    self.relation_dict[term.name][relation]["id"].append(value.id[0])
            if 'hasExactSynonym' in term.other:
                for elem in term.other['hasExactSynonym']:
                    # remove some words like cell or cells
                    querywords = elem.split()
                    resultwords = [word for word in querywords if word.lower() not in stopwords]
                    elem = ' '.join(resultwords)
                    elem = elem.lower()
                    self.relation_dict[term.name]['hasExactSynonym'].add(elem)
            if 'hasRelatedSynonym' in term.other:
                for elem in term.other['hasRelatedSynonym']:
                    querywords = elem.split()
                    resultwords = [word for word in querywords if word.lower() not in stopwords]
                    elem = ' '.join(resultwords)
                    elem = elem.lower()
                    self.relation_dict[term.name]['hasRelatedSynonym'].add(elem)
            if 'xref' in term.other:
                for elem in term.other['xref']:
                    self.relation_dict[term.name]['xref'].add(elem)
            if 'onProperty' in term.other and 'someValuesFrom' in term.other:
                clean_property_array = []
                for idx, elem in enumerate(term.other['onProperty']):
                    if '#' in term.other['onProperty'][idx] and "has_origin" not in term.other['onProperty'][idx]:
                        relation = str(elem).split("#")[1]
                        clean_property_array.append(relation)

                for idx, elem in enumerate(clean_property_array):
                    value_id_from = str(term.other['someValuesFrom'][idx]).split("_")[1]
                    value_name_from_array = str(term.other['someValuesFrom'][idx]).split("_")[0].split("/")
                    value_name_from = value_name_from_array[len(value_name_from_array) - 1]

                    querywords = value_name_from.split()
                    resultwords = [word for word in querywords if word.lower() not in stopwords]
                    value_name_from = ' '.join(resultwords)
                    value_name_from = value_name_from.lower()
                    # print(elem, value_name_from + ":" + value_id_from)
                    self.relation_dict[term.name][elem].append(value_name_from + ":" + value_id_from)
            # print(term.other, term.id, term)
            # print("\n")

        self.concate_part_of_is_a_relations(50)

    def get_nth_generation_relation(self, level, relation, dest, child, parent):
        """
        :param level: recursion depth is the depth of the graph
        :param relation: i.e. "is_a"
        :param child: given by the caller
        :param parent: the node were the childs should be added
        :param dest: position/key of the parent node were childs should be added
        """
        # remove some words like cell or cells
        stopwords = []
        querywords = child.split()
        resultwords = [word for word in querywords if word.lower() not in stopwords]
        child = ' '.join(resultwords)
        child = child.lower()
        if level == 0:
            return
        if "name" in self.relation_dict[child][relation]:
            for nthChild in self.relation_dict[child][relation]["name"]:
                self.relation_dict[parent][dest].add(nthChild)
                self.get_nth_generation_relation(level - 1, relation, dest, nthChild, parent)

    def concate_part_of_is_a_relations(self, level):
        """
        :param level: recursion depth is the depth of the graph
        concate information about is_a of n levels of the graph
        """
        for tissue, tissue_value in self.relation_dict.items():
            if "name" in tissue_value["is_a"]:
                for elem in tissue_value["is_a"]["name"]:
                    self.relation_dict[tissue]["is_a_childs"].add(elem)
                    self.get_nth_generation_relation(level, "is_a", "is_a_childs", elem, tissue)
            #if "name" in tissue_value["can_be"]:
            #    for elem in tissue_value["can_be"]["name"]:
            #        self.relation_dict[tissue]["is_a_childs"].add(elem)
            #        self.get_nth_generation_relation(level, "can_be", "is_a_childs", elem, tissue)



    def get_matching_node(self, query):
        """
        :param query: search disease in ont
        :return: tupel, where first is the similarity score and second pos. is the dict entry
        """
        cell_query = query.lower()
        clean_string = []
        for word in cell_query.split(" "):

            e = ''.join(e for e in word if e.isalnum() or e != "-")
            if e != "":
                clean_string.append(e)
        cell_query = " ".join(clean_string)
        best_match = [0, "", ""]
        tissue_in_ont = ""
        if cell_query in self.relation_dict:
            return [100, cell_query, ""]
        else:
            for disease, value in self.relation_dict.items():
                # if the first char does not match continue
                try:
                    if disease[0].lower() != cell_query[0].lower() or len(cell_query) < 3:
                        partial = 0
                    else:
                        partial = fuzz.ratio(disease, cell_query)
                except BaseException:
                    partial = 0
                for syn in value["hasRelatedSynonym"]:
                    try:
                        if syn[0].lower() != cell_query[0].lower() or len(cell_query) < 3:
                            partial_related_syn = 0
                        else:
                            partial_related_syn = fuzz.ratio(cell_query, syn)
                    except BaseException:
                        partial_related_syn = 0
                    # print("FuzzyWuzzy Ratio: ", fuzz.ratio(tissue, text), tissue)
                    # print("FuzzyWuzzy Ratio_PARTIAL: ", partial, tissue)
                    if best_match[0] < partial:
                        if disease.lower()[0] == cell_query[0] and len(cell_query) * 2 > len(disease):
                            best_match = [partial, disease, ""]
                    if best_match[0] < partial_related_syn:
                        if syn.lower()[0] == cell_query[0]:
                            best_match = [partial_related_syn, disease, syn]
                for syn in value["hasExactSynonym"]:
                    try:
                        if syn[0].lower() != cell_query[0].lower() or len(cell_query) < 3:
                            partial_related_syn = 0
                        else:
                            partial_related_syn = fuzz.ratio(cell_query, syn)
                    except BaseException:
                        partial_related_syn = 0
                    # print("FuzzyWuzzy Ratio: ", fuzz.ratio(tissue, text), tissue)
                    # print("FuzzyWuzzy Ratio_PARTIAL: ", partial, tissue)
                    if best_match[0] < partial:
                        if disease.lower()[0] == cell_query[0] and len(cell_query) * 2 > len(disease):
                            best_match = [partial, disease, ""]
                    if best_match[0] < partial_related_syn:
                        if syn.lower()[0] == cell_query[0]:
                            best_match = [partial_related_syn, disease, syn]
                if best_match[0] < partial:
                    if disease.lower()[0] == cell_query[0]:
                        best_match = [partial, disease, ""]
        return best_match

    def get_matching_node_token_set_ratio(self, query, threshold):
        """
        :param query: search term in ont
        :return: list of tupel, where first is the similarity score and second pos. is the dict entry
        Instead of just tokenizing the strings, sorting and then pasting the tokens back together,
        token_set_ratio performs a set operation that takes out the common tokens (the intersection) and
        then makes fuzz.ratio() pairwise comparisons between the following new strings:
        :param threshold: threshold
        """
        results = [[0, "", ""]]
        if query == "":
            return results
        # normalize string query
        search_query = query.lower()
        clean_string = []
        for word in search_query.split(" "):
            e = ''.join(e for e in word if e.isalnum() or e != "-")
            if e != "":
                clean_string.append(e)
        search_query = " ".join(clean_string)
        wnl = WordNetLemmatizer()
        tokens = [token.lower() for token in word_tokenize(search_query)]
        lemmatized_words = [wnl.lemmatize(token) for token in tokens]
        search_query = " ".join(lemmatized_words)
        for ont_term, value in self.relation_dict.items():
            if len(ont_term.split(" ")) > len(search_query.split(" ")) or len(ont_term) > len(search_query):
                partial = 0
            else:
                partial = fuzz.token_set_ratio(search_query, ont_term)

            for syn in value["hasRelatedSynonym"]:
                if len(syn.split(" ")) > len(search_query.split(" ")) or len(syn) > len(search_query):
                    continue
                partial_related_syn = fuzz.token_set_ratio(search_query, syn)
                if partial_related_syn > partial and partial_related_syn >= threshold:
                    results.append([partial_related_syn, ont_term, syn])

            for syn in value["hasExactSynonym"]:
                if len(syn.split(" ")) > len(search_query.split(" ")) or len(syn) > len(search_query):
                    continue
                partial_related_syn = fuzz.token_set_ratio(search_query, syn)
                if partial_related_syn > partial and partial_related_syn >= threshold:
                    results.append([partial_related_syn, ont_term, syn])

            if partial >= threshold:
                results.append([partial, ont_term, ""])
        return results

    @staticmethod
    def get_matching_node_st(query, relation_dict):
        """
        :param query: search disease in ont
        :return: tupel, where first is the similarity score and second pos. is the dict entry
        """
        cell_query = query.lower()
        clean_string = []
        for word in cell_query.split(" "):
            e = ''.join(e for e in word if e.isalnum() or e != "-")
            if e != "":
                clean_string.append(e)
        cell_query = " ".join(clean_string)
        wnl = WordNetLemmatizer()
        tokens = [token.lower() for token in word_tokenize(cell_query)]
        lemmatized_words = [wnl.lemmatize(token) for token in tokens]
        cell_query = " ".join(lemmatized_words)

        best_match = [0, "", ""]
        tissue_in_ont = ""
        if cell_query in relation_dict:
            return [100, cell_query, ""]
        else:
            for disease, value in relation_dict.items():
                # if the first char does not match continue
                try:
                    if disease[0].lower() != cell_query[0].lower() or len(cell_query) < 3:
                        partial = 0
                    else:
                        partial = fuzz.ratio(disease, cell_query)
                except BaseException:
                    partial = 0
                for syn in value["hasRelatedSynonym"]:
                    try:
                        if syn[0].lower() != cell_query[0].lower() or len(cell_query) < 3:
                            partial_related_syn = 0
                        else:
                            partial_related_syn = fuzz.ratio(cell_query, syn)
                    except BaseException:
                        partial_related_syn = 0
                    # print("FuzzyWuzzy Ratio: ", fuzz.ratio(tissue, text), tissue)
                    # print("FuzzyWuzzy Ratio_PARTIAL: ", partial, tissue)
                    if best_match[0] < partial:
                        if disease.lower()[0] == cell_query[0] and len(cell_query) * 2 > len(disease):
                            best_match = [partial, disease, ""]
                    if best_match[0] < partial_related_syn:
                        if syn.lower()[0] == cell_query[0]:
                            best_match = [partial_related_syn, disease, syn]
                for syn in value["hasExactSynonym"]:
                    try:
                        if syn[0].lower() != cell_query[0].lower() or len(cell_query) < 3:
                            partial_related_syn = 0
                        else:
                            partial_related_syn = fuzz.ratio(cell_query, syn)
                    except BaseException:
                        partial_related_syn = 0
                    # print("FuzzyWuzzy Ratio: ", fuzz.ratio(tissue, text), tissue)
                    # print("FuzzyWuzzy Ratio_PARTIAL: ", partial, tissue)
                    if best_match[0] < partial:
                        if disease.lower()[0] == cell_query[0] and len(cell_query) * 2 > len(disease):
                            best_match = [partial, disease, ""]
                    if best_match[0] < partial_related_syn:
                        if syn.lower()[0] == cell_query[0]:
                            best_match = [partial_related_syn, disease, syn]
                if best_match[0] < partial:
                    if disease.lower()[0] == cell_query[0]:
                        best_match = [partial, disease, ""]
        return best_match

    @staticmethod
    def get_matching_node_token_set_ratio_st(query, threshold, relation_dict):
        """
        :param query: search term in ont
        :return: list of tupel, where first is the similarity score and second pos. is the dict entry
        Instead of just tokenizing the strings, sorting and then pasting the tokens back together,
        token_set_ratio performs a set operation that takes out the common tokens (the intersection) and
        then makes fuzz.ratio() pairwise comparisons between the following new strings:
        :param threshold: threshold
        """
        results = [[0, "", ""]]
        if query == "":
            return results
        # normalize string query
        search_query = query.lower()
        clean_string = []
        for word in search_query.split(" "):
            e = ''.join(e for e in word if e.isalnum() or e != "-")
            if e != "":
                clean_string.append(e)
        search_query = " ".join(clean_string)
        wnl = WordNetLemmatizer()
        tokens = [token.lower() for token in word_tokenize(search_query)]
        lemmatized_words = [wnl.lemmatize(token) for token in tokens]
        search_query = " ".join(lemmatized_words)

        for ont_term, value in relation_dict.items():
            if len(ont_term.split(" ")) > len(search_query.split(" ")) or len(ont_term) > len(search_query):
                partial = 0
            else:
                partial = fuzz.token_set_ratio(search_query, ont_term)

            for syn in value["hasRelatedSynonym"]:
                if len(syn.split(" ")) > len(search_query.split(" ")) or len(syn) > len(search_query):
                    continue
                partial_related_syn = fuzz.token_set_ratio(search_query, syn)
                if partial_related_syn > partial and partial_related_syn >= threshold:
                    results.append([partial_related_syn, ont_term, syn])

            for syn in value["hasExactSynonym"]:
                if len(syn.split(" ")) > len(search_query.split(" ")) or len(syn) > len(search_query):
                    continue
                partial_related_syn = fuzz.token_set_ratio(search_query, syn)
                if partial_related_syn > partial and partial_related_syn >= threshold:
                    results.append([partial_related_syn, ont_term, syn])

            if partial >= threshold:
                results.append([partial, ont_term, ""])
        return results

class Phylogenetic_ontology(Disease_ontology):
     def __init__(self, location="ncbitaxon.owl"):
        self.ont = pronto.Ontology(location)
        self.relation_dict = dict()
        self.build_graph_doid()


if __name__ == '__main__':
    include_attr ={}
    include_attr["Host_DISEASE"] = Disease_ontology(location="2021_07_09__doid.obo")