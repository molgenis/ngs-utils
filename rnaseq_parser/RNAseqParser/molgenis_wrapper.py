###################################################################
#
# Molgenis python api client. 
#
####################################################################

import json
import requests
import os.path
from RNAseqParser import security
import time
import logging
from datetime import datetime
import timeit
from RNAseqParser import molgenis
class Connect_Molgenis():
    """This class only has __enter__ and __exit__ function to force use of with statement. This way the passwords saved to file can be cleaned up"""
    def __init__(self, server_url, remove_pass_file = True, new_pass_file = True, password_location = '~',log_file = 'molgenis.log', logging_level='DEBUG', logfile_mode = 'w', profile=True):
        """Pass all args and kwargs to the actual Connection class in __enter__()"""
        self.server_url = server_url
        self.remove_pass_file = remove_pass_file
        self.new_pass_file = new_pass_file
        self.password_location = password_location
        self.log_file = log_file
        self.logging_level = logging_level
        self.logfile_mode = logfile_mode
        self.profile = profile
    def __enter__(self):
        self.enter = True
        class Connection():
            """Actual Class with functionallity. Some simple methods for adding, updating and retrieving rows from Molgenis though the REST API
            BELOW IS OUTDATED DOCS
        
            Args:
                  server_url (string): The url to the molgenis server
                
            Example:
                # make a connection
                with open molgenis.Connect_Molgenis()('http://localhost:8080') as connection:
                    # add a row to the entity public_rnaseq_Individuals
                    connection.add_entity_row('public_rnaseq_Individuals',{'id':'John Doe','age':'26', 'gender':'Male'})
                    # get the rows from public_rnaseq_Individuals where gender = Male
                    print connection.get('public_rnaseq_Individuals',[{'field':'gender', 'operator':'EQUALS', 'value':'Male'}])['items'] 
                    # update row in public_rnaseqIndivduals where id=John Doe -> set gender to Female
                    connection.update_entity_row('public_rnaseq_Individuals',[{'field':'id', 'operator':'EQUALS', 'value':'John Doe'}], {'gender':'Female'})  
            """
        
            def __init__(self, server_url, remove_pass_file = True, new_pass_file = True, password_location = '~',log_file = 'molgenis.log', logging_level='ERROR', logfile_mode = 'w'):
                '''Initialize Python api to talk to Molgenis Rest API
                
                Args:
                    server_url (string):         The url to the molgenis server (ex: https://molgenis39.target.rug.nl/)
                    remove_pass_file (bool):     If True, remove the files containing the passwords after usage (def: True)
                    new_pass_file (str):         If file with password was not removed after last run, but still want to use a new password this run, set to True. Otherwise uses same password as last run (def: False)
                    password_location (string):  Folder where to put the password files in (def: ~)
                    log_file (string):           Path to write logfile with debug info etc to (def: molgenis.log)
                    logging_level (string):      The level of logging to use. See Python's `logging` manual for info on levels (def: DEBUG)
                    logfile_mode (string):       Mode of writing to logfile, e.g. w for overwrite or a for append, see `logging` manual for more details (def: w)
                '''
                # because errors in the __init__ function will not go to __exit__, make sure to clean up after error
                try:
                    # if no path is specified in the log_file name, it should be written in the same location where the script is called from,
                    # not from the location molgenis is located
                    if not os.sep in log_file:
                        log_file = os.getcwd()+os.sep+log_file
                    else:
                        # if there is a path in log_file, make sure that the folder exists
                        if not os.path.exists(os.path.dirname(log_file)):
                            raise OSError('Folder "'+str(os.path.dirname)+'" for writing the molgenis.log file does not exist, change log_file location')
                    logging.basicConfig(filename = log_file, filemode = logfile_mode)
                    logging.getLogger().addHandler(logging.StreamHandler())
                    self.logger = logging.getLogger(__name__)
                    self.logger.setLevel(level=getattr(logging, logging_level))
                    self.time_start = timeit.default_timer()
                    security.overwrite_passphrase_location(password_location)
                    if new_pass_file:
                        self.remove_pass_file = True
                        security.remove_secrets_file()
                    security.require_username('Username')
                    security.require_password('Password')
                    self.server_url = server_url
                    if not server_url.endswith('api/'):
                        if not server_url.endswith('/'):
                            server_url += '/'
                        server_url += 'api/'
                    self.session = molgenis.Session(server_url)
                    self.logger.debug('Trying to log in with data from '+str(security.PASSPHRASE_FILE) +' to: '+server_url+' with username: '+'*'*len(security.retrieve('Username'))+' password: '+'*'*len(security.retrieve('Password')))
                    try:
                        self.session.login(security.retrieve('Username'), security.retrieve('Password'))
                    except requests.exceptions.HTTPError as e:
                        if 'Not Found for url' in str(e):
                            self.logger.debug('login failed, trying again')
                            if server_url.startswith('http:'):
                                server_url = server_url.replace('http:','https:')
                                self.session = molgenis.Session(server_url)
                            elif server_url.startswith('https:'):
                                server_url = server_url.replace('https:','http:')
                                self.session = molgenis.Session(server_url)
                            self.logger.debug('Trying to log in with data from '+str(security.PASSPHRASE_FILE) +' to: '+server_url+' with username: '+'*'*len(security.retrieve('Username'))+' password: '+'*'*len(security.retrieve('Password')))
                            try:
                                self.session.login(security.retrieve('Username'), security.retrieve('Password'))
                            except requests.exceptions.HTTPError as e:
                                if 'Unauthorized for url' in str(e):
                                    raise requests.exceptions.HTTPError(str(e)+'\nInvalid username or password')
                                else:
                                    raise
                        elif 'Unauthorized for url' in str(e):
                            raise requests.exceptions.HTTPError(str(e)+'\nInvalid username or password')
                        else:
                            if len(security.retrieve('Username')) == 0 or len(security.retrieve('Password')) == 0:
                                raise requests.exceptions.HTTPError(str(e)+'\nError possibly because username or password is empty')
                            else:
                                raise
                    self.entity_meta_data = {}
                    self.column_meta_data = {}
                    self.added_rows = 0
                    self.added_files = 0
                    self.remove_pass_file = remove_pass_file
                except:
                    self.remove_password_files()
                    raise
                           
            def _sanitize_data(self, data, datetime_column, added_by_column):
                '''Remove None and empty lists from data. datetime_columnd = added_by or updated_by
                Args
                    ids_to_remove -- data with these IDs will be removed from the dataset
                '''
                data.update({datetime_column : str(datetime.now())})
                data.update({added_by_column : security.retrieve('Username')})
                # make all values str and remove if value is None or empty string
                data = {k: v for k, v in list(data.items()) if v!=None}
                data = dict([a, str(x)] for a, x in data.items() if len(str(x).strip())>0)
                return data
                        
            def _logging(self, entity, add_type):
                '''Add datetime and added by to entity row or file row
                
                entity (string): Name of the entity
                add_type (string): Either entity_row or file'''
                if add_type == 'entity_row': 
                    message = time.strftime('%H:%M:%S', time.gmtime(timeit.default_timer()-self.time_start))+' - Add row to entity '+entity+'. Total: '+str(self.added_rows)
                elif add_type == 'file':
                    message = time.strftime('%H:%M:%S', time.gmtime(timeit.default_timer()-self.time_start))+ ' - Add a file to '+entity+'. Total: '+str(self.added_files)
                else:
                    raise ValueError('add_type can only be entity_row or file')
                self.logger.debug(message)
                
            def _check_duplicate(self, entity, data):
                meta_data = self.get_entity_meta_data(entity)
                unique_attributes = []
                for attribute in meta_data['attributes']:
                    if meta_data['attributes'][attribute]['unique'] == True:
                        unique_attributes.append(attribute)
                id_list = []
                for attribute in unique_attributes:
                    if isinstance(data,dict):
                        if attribute in data:
                            id_list.append(data[attribute])
                    elif isinstance(data,(list,tuple)):
                        for data_ in data:
                            try:
                                id_list.append(data_[attribute])
                            except:
                                pass
                    else:
                        raise TypeError('data should be of type dict, list or tuple')
                duplicate_ids = []
                to_remove = []
                for result in self.session.get(entity,num=10000):
                    if attribute in result and result[attribute] in id_list:
                        if isinstance(data,dict):
                            return {},[result[self.get_id_attribute(entity)]]
                        duplicate_ids.append(result[attribute])
                        to_remove.append(id_list.index(result[attribute]))
                for index in sorted(to_remove, reverse=True):
                    del(data[index])   
                return data, duplicate_ids
                                         
            def add(self, entity, data, files = {}):
                '''Add one or multiple rows to an entity, or add a file
                
                Args:
                    entity (string): Name of the entity where row should be added
                    data (list or dict): To add multiple rows: list of dicts with Key = column name, value = column value
                                         To add one row: dict with Key = column name, value = column value
                    files -- dictionary containing file attribute values for the entity row.
                            The dictionary should for each file attribute map the attribute name to a tuple containing the file name and an input stream.

                Returns:
                    added_ids (list): List of IDs of the rows that got added
                    
                Example:
                    >>> with Connect_Molgenis('https://localhoost:8080') as connection:
                    >>>     print (connection.add('EntityName',[{'column_A':'row 1','column_B', 'row 1'},{'column_A':'row 2','column_B':'row 2'}])
                    >>> with Connect_Molgenis('https://localhoost:8080') as connection:
                    >>>>     print (connection.add('EntityName',{'column_A':'row 1','column_B', 'row 1'})
                    AAAACUGUI6T5KJXRMQK476QAAE
                    
                '''
                if not isinstance(data,list) and not isinstance(data,dict):
                    raise TypeError('data_list should be of type list or dict')
                if not isinstance(data,list):
                    data = [data]
                if len(data) == 0:
                    raise ValueError('data_list is an empty list, needs to contain a dictionary')
                #data,duplicate_ids = self._check_duplicate(entity, data)
                sanitized_data_list = [self._sanitize_data(data,'datetime_added','added_by') for data in data]
                # post to the entity with the json data
                if len(sanitized_data_list) == 1 and len(sanitized_data_list[0]) == 0:
                    return None
                try:
                    added_ids = self.session.add_all(entity, sanitized_data_list)
                except requests.exceptions.HTTPError as e:
                    try:
                        # if there is json() information, add it to the error message, as this gives the user information on actual problem
                        self.logger.debug(str(e.response.json()))
                    except AttributeError:
                        pass
                    raise
                self.added_rows += len(sanitized_data_list)
                #added_ids += duplicate_ids
                self._logging(entity,'entity_row')
                return added_ids

            def add_file(self, file_path, description, entity, data={}, file_name=None):
                '''Add a file to entity File. (io stream is possible but not implemented, see molgenis_api.py in archive for example)
                
                Args:
                    file_path (string): Path to the file to be uploaded
                    description (description): Description of the file
                    entity (string): Name of the entity to add the files to
                    data (dict): If extra columns have to be added, provide a dict with key column name, value value (def: None)
                    file_name (string): Name of the file. If None is set to basename of filepath (def: None)
                Returns:
                    file_id (string): ID if the file that got uploaded (for xref)
                    
                '''
                self._logging(entity,'file')
                data, added_id = self._check_duplicate(entity,data)
                if added_id:
                    return added_id
                data.update({'description': description})
                data = self._sanitize_data(data,'datetime_added','added_by')
                if not file_name:
                    file_name = os.path.basename(file_path)
                if not os.path.isfile(file_path):
                    self.logger.error('File not found: '+str(file_path))
                    raise IOError('File not found: '+str(file_path))
                added_id = self.session.add(entity, data=data,
                                            files={'attachment':(os.path.basename(file_path), open(file_path,'rb'))})
                self.added_files += 1
                return added_id

                       
            def get(self, entity, query=None,attributes=None, num=100, start=0, sortColumn=None, sortOrder=None):
                '''Get row(s) from entity with a query
                
                Args:
                    entity (string): Name of the entity where get query should be run on
                    query (list): List of dictionaries with as keys:values -> [{'field':column name, 'operator':'EQUALS', 'value':value}]
                    attributes (?): Attributes to return
                    num (int): Number of results to return
                    start (int): Page to start returning results from
                    sortColumn (string): Column name to sort on
                    sortOrder (string): ORder to sort in
                Returns:
                    result (dict): json dictionary of retrieve data
                '''
                entity = entity.strip()
                if query:
                    if len(query) == 0:
                        self.logger.error('Can\'t search with empty query')
                        raise ValueError('Can\'t search with empty query')
                    try:
                        items = self.session.get(entity, q = query, attributes=attributes, num=num, start=start, sortColumn=sortColumn, sortOrder=sortOrder)
                    except requests.exceptions.HTTPError as e:
                        try:
                            self.logger.debug(e.response.json())
                        except json.decoder.JSONDecodeError:
                            pass
                        raise
                else:
                    items = self.session.get(entity,attributes=attributes, num=num, start=start, sortColumn=sortColumn, sortOrder=sortOrder)         
                return items
          
            def update_entity_rows(self, entity, data, row_id = None, query_list=None):
                '''Update an entity row, either by giving the attribute id name and the id for the row to update, or a query for which row to update
            
                Args:
                    entity (string): Name of the entity to update
                    data (dict):  Key = column name, value = column value
                    id_attribute: The id_attribute name for the entity which you want to update the row
                    row_id: The row id value (from id_attribute)
                    query_list (list): List of dictionaries which contain query to select the row to update (see documentation of get())  (def:None)
                    validate_json (bool): If True, check if the given data keys correspond with the column names of entity. (def: False)
                '''
                data = self._sanitize_data(data,'datetime_last_updated','updated_by')
                id_attribute = self.get_id_attribute(entity)
                if row_id:
                    if query_list:
                        logging.warn('Both row_id and query_list set, will use only row_id')
                    try:
                        server_response = self.session.update(entity, row_id, data)
                    except requests.exceptions.HTTPError as e:
                        self.logger.debug(e.response.json())
                        raise
                    self.logger.debug(time.strftime('%H:%M:%S', time.gmtime(timeit.default_timer()-self.time_start))+' - Updated value of entity '+entity)
                    return server_response
                elif query_list:
                    queries = []
                    for query in query_list:
                        queries.append(self._sanitize_data(query,'datetime_last_updated','updated_by'))
                    entity_data = self.get(entity, queries)
                    if len(entity_data) == 0:
                        self.logger.error('Query returned 0 results, no row to update.')
                        raise Exception('Query returned 0 results, no row to update.')
                    for entity_item in entity_data:
                        server_response = self.session.update(entity, str(entity_item[str(id_attribute)]), data)
                        server_response = server_response
                        self.logger.debug(time.strftime('%H:%M:%S', time.gmtime(timeit.default_timer()-self.time_start))+' - Updated value of entity '+entity)
                    return server_response
                else:
                    raise ValueError('update_entity_rows function called without setting either row_id or query_list (one of the two needed to know which row to update)')
        
            def get_entity_meta_data(self, entity):
                '''Get metadata from entity
                
                Args:
                    entity (string): Name of the entity to get meta data of
                
                Returns:
                    result (dict): json dictionary of retrieve data
                '''
                if entity in self.entity_meta_data:
                    return self.entity_meta_data[entity]
                try:
                    entity_meta_data = self.session.get_entity_meta_data(entity)
                except requests.exceptions.HTTPError as e:
                    raise requests.exceptions.HTTPError(str(e)+\
                            '\nCheck if the package name in Config is correct and that an EMX is uploaded to Molgenis server '+str(self.server_url ))
                self.entity_meta_data[entity] = entity_meta_data
                return entity_meta_data
        
            def get_column_names(self, entity):
                '''Get the column names from the entity
                
                Args:
                    entity (string): Name of the entity to get column names of
                Returns:
                    meta_data(list): List with all the column names of entity
                '''
                entity_meta_data = self.get_entity_meta_data(entity)
                attributes = entity_meta_data['attributes']
                return list(attributes.keys())
            
            def get_id_attribute(self, entity):
                '''Get the id attribute name'''
                entity_meta_data = self.get_entity_meta_data(entity)
                return entity_meta_data['idAttribute']
            
            def get_column_meta_data(self, entity, column_name):
                '''Get the meta data for column_name of entity
                
                Args:
                    entity (string): Name of the entity 
                    column_name (string): Name of the column
                Returns:
                    List with all the column names of entity
                '''
                if entity+column_name in self.column_meta_data:
                    return self.column_meta_data[entity+column_name]
                server_response = self.session.get(self.api_v1_url+'/'+entity+'/meta/'+column_name)
                column_meta_data = server_response.json()
                self.column_meta_data[entity+column_name] = column_meta_data
                return column_meta_data
            
            def get_column_type(self, entity, column_name):
                column_meta_data = self.get_column_meta_data(entity, column_name)
                return column_meta_data['fieldType']
                   
            def delete_all_rows_of_all_entities(self, package):
                '''Delete all entities of package
                
                Args:
                    package (string): Package for which to delete all entities. (def: None)
                '''
                if not package:
                    self.logger.error('package can\'t be None, is '+str(package))
                    raise AttributeError('package can\'t be None, is '+str(package))
                server_response = self.get_all_entity_data()
                for entity in server_response.json()['items']:
                    entity = entity['fullName']
                    if package in entity and not bool(entity['abstract']):
                        self.logger.info('Deleting all rows from',entity)
                        try:
                            self.delete_all_entity_rows(entity)
                        except Exception as e:
                            self.logger.warning(str(e))
            
            def delete_all_entity_rows(self,entity):
                '''delete all entity rows'''
                items = self.get(entity,num=10000)
                while len(items) > 0:
                    self.delete_entity_data(entity,items)
                    items = self.get(entity,num=10000)
            def delete_entity_rows(self, entity, query):
                '''delete entity rows
            
                Args:
                    entity (string): Name of the entity to update
                    query (list): List of dictionaries which contain query to select the row to update (see documentation of get())
                '''
                items = self.get(entity, query)
                if len(items) == 0:
                    self.logger.error('Query returned 0 results, no row to delete.')
                    raise Exception('Query returned 0 results, no row to delete.')
        
            def delete_entity_data(self, entity, items):
                '''delete entity data
                
                Args:
                    entity_data (dict): A dictionary with at least key:"items", value:<dict with column IDs>. All items in this dict will be deleted
                    entity (string): Name of entity to delete from
                    query_used (string): Incase entity_data was made with a query statement, the query used can be given for more detailed error prints (def: None)
                '''
                id_attribute = self.get_id_attribute(entity)
                for rows in items:
                    row_id = rows[id_attribute]
                    try:
                        self.logger.debug('Deleted row from: '+str(entity))
                        self.session.delete(entity.strip(),row_id)
                    except requests.exceptions.HTTPError as e:
                        self.logger.debug(e.response.json())
                        raise
            
            def remove_password_files(self):
                if self.remove_pass_file:
                    security.remove_secrets_file()
                
        self.molgenis_connection_obj = Connection(self.server_url,
                                                  remove_pass_file = self.remove_pass_file,
                                                  new_pass_file = self.new_pass_file,
                                                  password_location = self.password_location,
                                                  log_file = self.log_file,
                                                  logging_level = self.logging_level,
                                                  logfile_mode = self.logfile_mode)
        return self.molgenis_connection_obj
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.molgenis_connection_obj.remove_password_files()

