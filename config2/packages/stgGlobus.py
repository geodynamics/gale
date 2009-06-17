import os
from Package import Package

class stgGlobus(Package):

    def gen_locations(self):
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib'])
        yield ('/usr/local', ['/usr/local/include/libxml2'], ['/usr/local/lib'])

    def gen_envs(self, loc):
        libs = "globus_io_gcc32dbg globus_ftp_control_gcc32dbg globus_ftp_client_gcc32dbg globus_gass_transfer_gcc32dbg globus_wsrf_core_tools_gcc32dbg globus_c_wsrf_resource_gcc32dbg globus_wsa_bindings_gcc32dbg globus_ds_bindings_gcc32dbg globus_wsu_bindings_gcc32dbg globus_wsp_bindings_gcc32dbg globus_wst_bindings_gcc32dbg globus_wsse_bindings_gcc32dbg globus_wsseu_bindings_gcc32dbg globus_wsc_bindings_gcc32dbg globus_wsbf_bindings_gcc32dbg globus_wsnt_bindings_gcc32dbg globus_wsrl_bindings_gcc32dbg globus_wsrp_bindings_gcc32dbg globus_wssg_bindings_gcc32dbg globus_notification_consumer_bindings_gcc32dbg globus_subscription_manager_bindings_gcc32dbg globus_ws_messaging_gcc32dbg globus_usage_gcc32dbg globus_xio_gcc32dbg gssapi_error_gcc32dbg globus_gss_assist_gcc32dbg globus_gssapi_gsi_gcc32dbg globus_gsi_proxy_core_gcc32dbg globus_gsi_credential_gcc32dbg globus_gsi_callback_gcc32dbg globus_oldgaa_gcc32dbg globus_gsi_sysconfig_gcc32dbg globus_gsi_cert_utils_gcc32dbg globus_openssl_gcc32dbg globus_openssl_error_gcc32dbg globus_callout_gcc32dbg globus_proxy_ssl_gcc32dbg globus_common_gcc32dbg ssl_gcc32dbg crypto_gcc32dbg xml2_gcc32dbg ltdl_gcc32dbg"
        for env in Package.gen_envs(self, loc):
            env['pkg_headers'] = [os.path.join('libxml', 'parser.h')]
            if self.find_libraries(loc[2], 'libxml'):
                env.PrependUnique(LIBS=libs.split())
                yield env
