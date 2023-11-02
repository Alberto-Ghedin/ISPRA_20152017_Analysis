import pymysql
from itertools import chain
import pandas as pd

class Database(): 

    def __init__(self, con, name, structure=None):
        self.connection = con
        self.cur = self.connection.cursor()
        self.name = name
        #structure is ditionary {table_name : columns}
        self.structure = structure

    def create_cur(self): 
        if self.cur is None: 
            self.cur = self.connection.cursor()


    def switch_to_db(self): 
        self.create_cur()
        self.cur.execute(f"USE {self.name};")


    def create_table(self, dict_table_colums):
        self.create_cur()
        for table_name, columns in dict_table_colums.items(): 
            # STRING CONTAINING COLUMN INFORMATION
            table_content = "{} {}, "*len(columns)
            # REMOVE LAST BLANK SPACE AND COMMA
            table_content = f"({table_content[:-2]});"
            self.cur.execute(f"CREATE TABLE IF NOT EXISTS {table_name}" + table_content.format(*chain(*columns.items())))
        
    
    def create_db(self):
        self.create_cur()
        self.cur.execute(f"CREATE DATABASE IF NOT EXISTS {self.name};")
        
        self.switch_to_db()

        # CREATING TABLES
        if self.structure is None: 
            raise ValueError("Structure must be specified!")
        self.create_table(self.structure)


    def rename_table(self, dict_new_old_table_name=None): 
        #dict_new_old_table_name = {old_name : new_name}
        self.switch_to_db()
        self.create_cur()
        if dict_new_old_table_name: 
            new_name_query = "RENAME TABLE {} TO {};"*len(dict_new_old_table_name.keys())
            self.cur.execute(new_name_query.format(*chain(*dict_new_old_table_name.items())))


    def alter_table(self, dict_old_new_col_name=None, add_cols=None, modify_cols=None, modify_key=None): 
        #dict_new_old_table_name = {old_name : new_name}
        #add_cols = {table_name : {cols_name : properties, ... }, ...}
        #add_table = {table_name : {cols_name : properties, ... }, ...}
        #modify_cols = {table_name : {cols_name : properties, ... }, ...}
        #modify_key = {table_name : {cols_name : properties, ... }, ...}

        if dict_old_new_col_name: 
            new_name_query = ""
            for table_name, old_new_cols in dict_old_new_col_name: 
                new_name_query = f"ALTER TABLE {table_name}" + " RENAME COLUMN {} TO {}"*len(old_new_cols)
                new_name_query = (new_name_query[:-1] + ";").format(*chain(*old_new_cols.items()))
            self.cur.execute(new_name_query)
        if add_cols: 
            add_cols_query = ""
            for table_name, columns in add_cols.items():
                add_cols_query += f"ALTER TABLE {table_name}" + " ADD COLUMN {} {},"*len(columns)
                add_cols_query = (add_cols_query[:-1] + ";").format(*chain(*columns.items()))
            self.cur.execute(add_cols_query)
        if modify_cols: 
            mod_cols_query = ""
            for table_name, columns in modify_cols.items(): 
                mod_cols_query += f"ALTER TABLE {table_name}" + " MODIFY COLUMN {} {},"*len(columns)
                mod_cols_query = (mod_cols_query[:-1] + ";").format(*chain(*columns.items()))
            self.cur.execute(mod_cols_query)
        if modify_key: 
            mod_key_query = ""
            for table_name, columns in modify_key.items(): 
                for column_name, new_property in columns.items(): 
                    self.cur.execute(
                        (
                            "SELECT CONSTRAINT_NAME, REFERENCED_TABLE_NAME, REFERENCED_COLUMN_NAME " 
                            "FROM INFORMATION_SCHEMA.KEY_COLUMN_USAGE "
                            "WHERE COLUMN_NAME = %s AND "
                            "TABLE_SCHEMA = %s AND "
                            "REFERENCED_TABLE_NAME IS NOT NULL"
                            ),
                        (column_name, self.name)
                    )

                    costraint_name, ref_table_name, ref_col_name = self.cur.fetchone()
                    self.cur.execute(f"ALTER TABLE {table_name} DROP FOREIGN KEY {costraint_name}")
                    self.cur.execute(
                        f"ALTER TABLE {table_name} ADD CONSTRAINT FOREIGN KEY"
                        + f"({column_name}) REFERENCES {ref_table_name}({ref_col_name})" 
                        + f"{new_property}"
                     )

    def __create_join_part(
            self, 
            join_type_table_names : dict[list[str]], 
            join_type_select : dict[list[list[str]]],  
            join_type_on : dict[list[dict]], 
            query_select
    ) -> tuple[str]: 
        query_join = ""

        for (type, join_table_names), join_select, join_on in zip(join_type_table_names.items(), join_type_select.values(), join_type_on.values()): 
            for indx, name in enumerate(join_table_names):
                        for col in join_select[indx]:  
                            if not col == "": 
                                query_select  += f", {name}.{col}"
                        query_join += " {} JOIN {} ON {} = {}".format(type, name, *join_on[indx].keys(), *join_on[indx].values())
            
        return query_join, query_select


    def extract_data(self,
                    main_table : str, 
                    selected_cols : list[str],  
                    where_clauses : dict = None, 
                    join_type_table_names : dict[list[str]] = None, 
                    join_type_select : dict[list[list[str]]] = None,  
                    join_type_on : dict[list[dict]] = None, 
                    print_query=False, 
                    distinct=False) -> pd.DataFrame:
        """
            function to extract data from sql database.

            where_clauses dictionary must be of the form 
            {
                "=" : [<colum_name_1>, <value_1>, <colum_name_2>, <value_2>, ...],
                *analogous for > < *, 
                "BETWEEN" : [<colum_name_1>, <lower_bound_1>, <upper_bound_1>, <colum_name_2>, <lower_bound_2>, <upper_bound_2>, ...]
            }

            if in the join not column is requested, than a "" should be put in the corresponding position of the join_select list
        """
        with self.connection.cursor(pymysql.cursors.DictCursor) as cur:
            query_select = "SELECT" 
            if distinct: 
                query_select += " DISTINCT"
            for item in selected_cols[:-1]:
                query_select += f" {main_table}.{item},"
            query_select += f" {main_table}.{selected_cols[-1]}"
            
            query_join = ""
            if join_type_table_names:
                query_join, query_select = self.__create_join_part(join_type_table_names, join_type_select, join_type_on, query_select)

            query_select += f" FROM {main_table}"

            if where_clauses: 
                query_where = " WHERE "
                if "BETWEEN" in  where_clauses.keys(): 
                    pairs = where_clauses["BETWEEN"]
                    N = int((len(pairs) / 3))
                    query_where += ("{} {} {} AND {} AND " * N).format(
                                                                *chain(
                                                                        *zip(pairs[0::3], ["BETWEEN"] * N, pairs[1::3], pairs[2::3])
                                                                        )
                                                                )
                for operator, pairs in filter(lambda x : x[0].upper() != "BETWEEN", where_clauses.items()): 
                    N = int((len(pairs) / 2))
                    query_where += ("{} {} {} AND " * N).format(
                                                                *chain(
                                                                        *zip(pairs[0::2], [operator] * N, pairs[1::2])
                                                                        )
                                                                )
                query_where = "".join(query_where.rsplit(" AND ", 1))
            else: 
                query_where = ""
            
            if print_query:
                print(cur.mogrify(query_select + query_join + query_where))
            cur.execute(query_select + query_join + query_where)
        return pd.DataFrame(cur.fetchall())
            

    def where_query(self, dict_target_col, print_query=False):
        self.create_cur()
        target_value = {}
        for new_col_name, (table, col, dict_values) in list(dict_target_col.items()):
            query_where = f"SELECT {col} AS {new_col_name} FROM {table} WHERE "
            # considering multiple keys for the selection
            keys = list(dict_values.keys())
            for key in keys[:-1]:
                query_where += f"{key}=%s AND "
            query_where += f"{keys[-1]}=%s;"

            # copying query for all the data
            query_where = query_where * \
                len(dict_values[keys[0]])
            if print_query: 
                print(self.cur.mogrify(query_where, 
                    tuple(
                        chain(*zip(*dict_values.values()))
                        )
                    )
                )
            self.cur.execute(query_where, 
                    tuple(
                        chain(*zip(*dict_values.values()))
                        )
                    )
            try: 
                target_value[new_col_name] = list(self.cur.fetchone())
            except: 
                if self.cur.fetchone() is None: 
                    target_value[new_col_name] = [self.cur.fetchone()]
            
            while self.cur.nextset() == True:
                try: 
                    target_value[new_col_name] += list(self.cur.fetchone())
                except: 
                    if self.cur.fetchone() is None: 
                        target_value[new_col_name] += [self.cur.fetchone()]
        return target_value


    def insert_into_table(self, Table_name, dict_cols_values=None, dict_forekeys=None, execute_unique=True, print_query=False, print_where_query=False):
        query_insert = f"INSERT INTO {Table_name} ("
        query_select = " SELECT "
        self.create_cur()
        self.switch_to_db()
        
        if dict_forekeys:
            ref_value = self.where_query(dict_forekeys, print_query=print_where_query)
            for key, value in ref_value.items():
                if len(value) == 1:
                    ref_value[key] = value * \
                        len(list(dict_cols_values.values())[0])
            try: 
                dict_cols_values.update(ref_value)
            except: 
                if dict_cols_values == None: 
                    dict_cols_values = ref_value

        # Second, specify insert and select entrance
        for col in dict_cols_values.keys():
            query_insert += col + ","
            query_select += "%s,"
        query_insert = query_insert[:-1] + ")"
        query_select = query_select[:-1]

        # where clause
        # find unique
        with self.connection.cursor(pymysql.cursors.DictCursor) as cur:
            if execute_unique:
                query_unique = "SELECT DISTINCT COLUMN_NAME FROM INFORMATION_SCHEMA.KEY_COLUMN_USAGE AS t1 INNER JOIN" \
                    + " INFORMATION_SCHEMA.TABLE_CONSTRAINTS AS t2 ON t1.CONSTRAINT_NAME= t2.CONSTRAINT_NAME WHERE CONSTRAINT_TYPE='UNIQUE'" \
                    + f" AND t2.TABLE_NAME='{Table_name}' AND t1.TABLE_NAME = '{Table_name}' AND t1.TABLE_SCHEMA='{self.name}'"
                cur.execute(query_unique)
                result = cur.fetchall()
                col_uniques = []
                for col in result:
                    col_uniques += list(col.values())
                query_exists = f" WHERE NOT EXISTS (SELECT * FROM {Table_name} WHERE "
                
                for unique in col_uniques[:-1]:
                    query_exists += unique + "=%s AND "
                query_exists += col_uniques[-1] + "= %s LIMIT 1);"
                
                query_final = query_insert + query_select + query_exists
                
                for i, datas in enumerate(zip(*list(dict_cols_values.values()))):
                    datas = list(datas)
                    for unique in col_uniques:
                        datas.append(dict_cols_values[unique][i]) 
                    if print_query:
                        print(cur.mogrify(query_final, tuple(datas)))
                    cur.execute(query_final, tuple(datas))
            else:
                query_final = query_insert + query_select
                for i, datas in enumerate(zip(*list(dict_cols_values.values()))):
                    datas = list(datas)
                    if print_query:
                        print(cur.mogrify(query_final, tuple(datas)))
                    cur.execute(query_final, tuple(datas))


    def Show_constraints(self, all_const : bool=False, where_clauses : dict = None) -> pd.DataFrame: 
        #part of the query common to all cases
        initial_part = ("SELECT "
                        "TABLE_NAME," 
                        "COLUMN_NAME,"
                        "CONSTRAINT_NAME,"
                        "REFERENCED_TABLE_NAME,"
                        "REFERENCED_COLUMN_NAME "
                        "FROM "
                            "INFORMATION_SCHEMA.KEY_COLUMN_USAGE "
                        "WHERE "
                            "TABLE_SCHEMA = %s")
        if all_const: 
            with self.connection.cursor(pymysql.cursors.DictCursor) as cur:
                cur.execute(initial_part, self.name)
                return pd.DataFrame(cur.fetchall())

        for constraint in where_clauses.keys(): 
            initial_part += f" AND {constraint}=%s"        
            
        with self.connection.cursor(pymysql.cursors.DictCursor) as cur:
            cur.execute(initial_part, (self.name, *where_clauses.values()))
            return pd.DataFrame(cur.fetchall())
        

    def delete_duplicate_from_table(self, table_name : str, print_query : bool = False): 
        self.cur.execute(f"SHOW COLUMNS FROM {table_name}")
        columns = self.cur.fetchall()
        column_names = list(zip(*columns))[0][1:]
        delete_query = f"DELETE t1 FROM {table_name} t1 INNER JOIN {table_name} t2 WHERE t1.id < t2.id AND "
        for name in column_names[:-1]:
            delete_query += f"t1.{name} = t2.{name} AND "
        delete_query += f"t1.{column_names[-1]} = t2.{column_names[-1]}"

        if print_query:
            print(self.cur.mogrify(delete_query))
        self.cur.execute(delete_query)

    
    def delete_duplicate_from_multiple_tables(self, table_names : list[str], print_query : bool = False): 
        for table_name in table_names: 
            self.delete_duplicate_from_table(table_name, print_query=print_query)
